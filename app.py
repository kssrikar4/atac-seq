import streamlit as st
import pandas as pd
import numpy as np
import glob
import os
import re
import plotly.express as px
import plotly.graph_objects as go

st.set_page_config(
    page_title="Chromatin Accessibility Analysis",
    page_icon="🧬",
    layout="wide"
)


st.markdown("""
    <style>
        .main-header { font-size: 2.2rem; font-weight: 700; margin-bottom: 0.5rem; }
        .sub-header { font-size: 1.1rem; color: #64748B; margin-bottom: 2rem; }
        .metric-card { background-color: #F8FAFC; border: 1px solid #E2E8F0; padding: 1rem; border-radius: 0.5rem; }
    </style>
""", unsafe_allow_html=True)

st.markdown('<div class="main-header">Chromatin Accessibility Analysis</div>', unsafe_allow_html=True)
st.markdown('<div class="sub-header">Interactive Epigenomics & Chromatin Structure Interpretation</div>', unsafe_allow_html=True)
st.markdown("---")

@st.cache_data
def parse_pipeline_data():
    frip_files = glob.glob("./qc/*_frip.txt")
    peak_files = glob.glob("./peaks/*/*_peaks.narrowPeak")
    samples = {}
    for f in frip_files:
        sample_id = os.path.basename(f).replace("_frip.txt", "")
        if sample_id not in samples:
            samples[sample_id] = {}
        try:
            with open(f, 'r') as file:
                samples[sample_id]['frip'] = float(file.read().strip())
        except Exception:
            samples[sample_id]['frip'] = 0.0
    for f in peak_files:
        sample_id = os.path.basename(os.path.dirname(f))
        if sample_id not in samples:
            samples[sample_id] = {}
        try:
            df = pd.read_csv(f, sep='\t', header=None, names=[
                'chrom', 'start', 'end', 'name', 'score', 'strand', 
                'signalValue', 'pValue', 'qValue', 'peak'
            ])
            samples[sample_id]['peaks_df'] = df
            samples[sample_id]['peak_count'] = len(df)
        except Exception:
            pass
    return samples

analysis_mode = st.sidebar.radio(
    "Select Operating Mode",
    ["Pipeline Directory Parse", "Direct Custom File Upload"]
)

selected_sample = None
sample_frip = None
peaks_df = None
total_peaks = 0

if analysis_mode == "Pipeline Directory Parse":
    pipeline_data = parse_pipeline_data()
    if not pipeline_data:
        st.warning("No pipeline-generated data detected in standard paths (./qc/ or ./peaks/). Switch to 'Direct Custom File Upload' or run the upstream pipeline first.")
        st.stop()
    
    selected_sample = st.sidebar.selectbox("Select Target Sample for Profiling", list(pipeline_data.keys()))
    sample_frip = pipeline_data[selected_sample].get('frip', 0.0)
    peaks_df = pipeline_data[selected_sample].get('peaks_df', None)
    total_peaks = pipeline_data[selected_sample].get('peak_count', 0)
    
else:
    st.sidebar.header("Direct File Upload")
    uploaded_peaks = st.sidebar.file_uploader(
        "Upload narrowPeak Files (.narrowPeak / .bed / .txt)", 
        type=["narrowPeak", "bed", "txt", "tsv"], 
        accept_multiple_files=True
    )
    
    selected_sample = "Uploaded Dataset"
    if uploaded_peaks:
        peak_list = []
        for i, f in enumerate(uploaded_peaks):
            try:
                df_temp = pd.read_csv(f, sep='\t', header=None, comment='#')
                if df_temp.shape[1] >= 3:
                    # Map standard BED columns, leave remaining columns generic
                    cols = ['chrom', 'start', 'end']
                    extra_cols = [f'col_{idx}' for idx in range(3, df_temp.shape[1])]
                    
                    # Map narrowPeak names if schema length matches
                    if df_temp.shape[1] == 10:
                        cols += ['name', 'score', 'strand', 'signalValue', 'pValue', 'qValue', 'peak']
                    else:
                        cols += extra_cols
                        
                    df_temp.columns = cols
                    df_temp['sample_origin'] = f.name
                    peak_list.append(df_temp)
            except Exception as e:
                st.sidebar.error(f"Error reading {f.name}: {str(e)}")
        if peak_list:
            peaks_df = pd.concat(peak_list, ignore_index=True)
            total_peaks = len(peaks_df)

if peaks_df is not None and 'length' not in peaks_df.columns:
    peaks_df['length'] = peaks_df['end'] - peaks_df['start']

col1, col2 = st.columns([1, 3])

with col1:
    st.markdown("### QC Metrics")
    if sample_frip is not None:
        st.metric("FRiP Score", f"{sample_frip:.4f}")
    st.metric("Active Peak Count", f"{total_peaks:,}")
    
    if peaks_df is not None:
        mean_len = peaks_df['length'].mean()
        median_len = peaks_df['length'].median()
        st.metric("Mean Peak Width", f"{mean_len:.1f} bp")
        st.metric("Median Peak Width", f"{median_len:.1f} bp")
    
    st.markdown("---")
    st.markdown("### Diagnostic Status")
    if sample_frip is not None:
        if sample_frip >= 0.20:
            st.success("**EXCELLENT**\nHighly enriched open chromatin with minimal background noise.")
        elif sample_frip >= 0.15:
            st.info("**ACCEPTABLE**\nStandard sample enrichment. Reliable downstream signals.")
        else:
            st.warning("**POOR ENRICHMENT**\nEnzymatic shearing or high cellular mortality. Inspect signal-to-noise ratio.")
    else:
        st.info("Direct File Import. Perform qualitative control assessments via tab panels.")

with col2:
    tab1, tab2, tab3, tab4 = st.tabs([
        "Chromatin Architecture Plots", 
        "Genomic Landscape Mapping", 
        "Signal & Significance Analysis",
        "Diagnostic Metrics Guide"
    ])
    
    with tab1:
        st.markdown("#### Physical Chromatin Fingerprint (Peak Width Distribution)")
        if peaks_df is not None:
            fig_hist = px.histogram(
                peaks_df, 
                x="length", 
                nbins=100, 
                color_discrete_sequence=["#3B82F6"],
                labels={'length': 'Peak Width (bp)'},
                title="Nucleosome-Free Regions (NFR) vs. Nucleosome Spacing"
            )
            fig_hist.add_vline(x=147, line_dash="dash", line_color="#EF4444", 
                             annotation_text="Nucleosomal Core (147bp)", annotation_position="top right")
            fig_hist.update_layout(
                plot_bgcolor="rgba(0,0,0,0)",
                paper_bgcolor="rgba(0,0,0,0)",
                margin=dict(l=20, r=20, t=40, b=20),
                xaxis=dict(gridcolor="#F1F5F9", range=[0, 1000]),
                yaxis=dict(gridcolor="#F1F5F9", title="Count")
            )
            st.plotly_chart(fig_hist, use_container_width=True)
            
            st.markdown("""
            **Structural Properties of Peaks:**
            * **< 147 bp (Left of line):** Nucleosome-Free Regions (NFR) typical of active promoters and accessible TF binding clusters.
            * **> 147 bp (Right of line):** Nucleosomal DNA fragments (mono, di, multi-nucleosomes) where linker regions were selectively targeted by Tn5.
            """)
        else:
            st.warning("Please load or import peak files to display architectural profiles.")

    with tab2:
        st.markdown("#### Chromosomal Peak Abundance")
        if peaks_df is not None:
            chrom_counts = peaks_df['chrom'].value_counts().reset_index()
            chrom_counts.columns = ['Chromosome', 'Peaks Count']
            
            chrom_order = sorted(chrom_counts['Chromosome'].unique(), key=lambda x: int(re.sub(r'\D', '', x)) if re.sub(r'\D', '', x).isdigit() else 99)
            chrom_counts['Chromosome'] = pd.Categorical(chrom_counts['Chromosome'], categories=chrom_order, ordered=True)
            chrom_counts = chrom_counts.sort_values('Chromosome')

            fig_chrom = px.bar(
                chrom_counts, 
                y="Chromosome", 
                x="Peaks Count", 
                orientation="h",
                color="Peaks Count",
                color_continuous_scale="Blues",
                title="Accessibility Distribution Across Chromosomes"
            )
            fig_chrom.update_layout(
                plot_bgcolor="rgba(0,0,0,0)",
                paper_bgcolor="rgba(0,0,0,0)",
                margin=dict(l=20, r=20, t=40, b=20),
                xaxis=dict(gridcolor="#F1F5F9"),
                yaxis=dict(autorange="reversed")
            )
            st.plotly_chart(fig_chrom, use_container_width=True)
        else:
            st.warning("Upload peak files to resolve landscape maps.")

    with tab3:
        st.markdown("#### Peak Signal Strength vs. Statistical Significance")
        if peaks_df is not None and 'signalValue' in peaks_df.columns and 'qValue' in peaks_df.columns:
            plot_df = peaks_df if len(peaks_df) <= 5000 else peaks_df.sample(n=5000, random_state=42)
            
            fig_scatter = px.scatter(
                plot_df, 
                x="signalValue", 
                y="qValue",
                color="length",
                color_continuous_scale="Viridis",
                labels={'signalValue': 'Signal Fold Enrichment', 'qValue': '-log10(q-value)'},
                title=f"Peak Significance Diagnostics ({len(plot_df):,} Coordinates Plotted)",
                hover_data=['chrom', 'start', 'end']
            )
            fig_scatter.update_layout(
                plot_bgcolor="rgba(0,0,0,0)",
                paper_bgcolor="rgba(0,0,0,0)",
                margin=dict(l=20, r=20, t=40, b=20),
                xaxis=dict(gridcolor="#F1F5F9"),
                yaxis=dict(gridcolor="#F1F5F9")
            )
            st.plotly_chart(fig_scatter, use_container_width=True)
            st.info("Top right elements represent high-enrichment, high-significance regions of active accessibility.")
        else:
            st.warning("Statistical plotting requires files containing enrichment variables (`signalValue` and `qValue`), e.g., standard .narrowPeak format.")

    with tab4:
        st.markdown("#### Diagnostic References Summary")
        
        col_ref1, col_ref2 = st.columns(2)
        with col_ref1:
            st.markdown("""
            ##### **Fraction of Reads in Peaks (FRiP)**
            $$\\text{FRiP} = \\frac{N_{\\text{reads in peaks}}}{N_{\\text{total mapped reads}}}$$
            * **High (>0.20):** Strong targeted transposase cuts in open, functional chromatin.
            * **Low (<0.15):** High background noise caused by non-specific genomic cuts or cell damage.
            """)
        with col_ref2:
            st.markdown("""
            ##### **Normalizations (CPM)**
            $$\\text{Normalized Score} = \\frac{\\text{Raw Bin Count} \\times 10^6}{N_{\\text{total mapped reads}}}$$
            * Normalization is standard across tracks to allow quantitative comparative signal analysis across divergent sample runs.
            """)
        
        st.markdown("---")
        st.markdown("""
        ##### **Downstream Applications Ready for Launch:**
        1. **Motif Enrichment & Transcription Factor Footprinting:** Identify localized target factor bind sites matching open footprints.
        2. **Differential Peaks Analyses (DESeq2/DiffBind):** Map differences in accessibility profiles across distinct operational conditions or time points.
        """)
