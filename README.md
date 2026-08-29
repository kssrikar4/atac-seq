# ATAC-seq

This guide is designed to take you from raw, unaligned sequencing reads to high-quality peaks and clean, normalized visualization tracks. In addition to the step-by-step commands, this repository includes two automated components designed to streamline your workflow and expand your analytics:

* **`run_pipeline.sh`**: A shell script that automates the entire end-to-end command-line processing (data download, QC, alignment, filtering, peak calling, BigWig generation, and quality reporting) with a single command.
* **`app.py`**: A interactive Streamlit dashboard. It acts as a diagnostic and biological report interface, allowing you to load your pipeline's outputs or directly upload independent `.bw` and `.narrowPeak` files for instant, data-driven structural and functional genomics insights.

## 1. Conda Environment Installation

To ensure absolute reproducibility and avoid polluting your global system paths, we create an isolated Conda environment.

```bash
# Create the environment inside the local directory './env'
conda create --prefix ./env -c bioconda -c conda-forge sra-tools samtools bowtie2 fastqc fastp macs3 deeptools bedtools subread pigz -y

# Activate the local environment
conda activate ./env

```

## 2. Directory Setup

```bash
mkdir -p ./refs ./raw ./fastq ./qc ./bam ./peaks ./bw ./tmp

```

### Directory Structure Map

* `./env`: Isolated Conda virtual environment containing all workflow dependencies.
* `./refs`: Reference genomes, indices, and genomic exclusion zones (blacklists).
* `./raw`: Uncompressed or compressed raw sequencing archive files.
* `./fastq`: De-multiplexed sequence data.
* `./qc`: Diagnostic reports (FastQC, fastp data) and processing metrics.
* `./bam`: Alignments, sorted, deduplicated, and clean sequence files.
* `./peaks`: Region coordinates marking open chromatin locations.
* `./bw`: Normalized signal tracks for genome browser viewing.
* `./tmp`: Large intermediate sorting cache files.


## 3. Setting Up the Genomic Foundations

To analyze where your reads belong, you need to map them back to a reference genome. This section sets up the genomic templates, indices for alignment acceleration, and an **ENCODE Blacklist** to discard known artifactual regions (e.g., highly repetitive regions that give false-positive peaks).

```bash
# 1. Download human reference genome (GRCh38)
wget -P ./refs/ https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GCA_000001405.15_GRCh38_genomic.fna.gz
gzip -d ./refs/GCA_000001405.15_GRCh38_genomic.fna.gz

# 2. Build the Bowtie2 index (This step can take 15–30 minutes)
bowtie2-build --threads 12 ./refs/GCA_000001405.15_GRCh38_genomic.fna ./refs/GRCh38_index

# 3. Index the fasta file and extract chromosome sizing profiles
samtools faidx ./refs/GCA_000001405.15_GRCh38_genomic.fna
cut -f1,2 ./refs/GCA_000001405.15_GRCh38_genomic.fna.fai > ./refs/GRCh38.chrom.sizes

# 4. Download the ENCODE Blacklist to filter out false signal peaks
wget -P ./refs/ https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz
gzip -d ./refs/ENCFF356LFX.bed.gz

```

## 4. Data Acquisition & Processing

Depending on your data source, you will handle your inputs in one of two ways. Choose the strategy below that matches your project setup.

### Mode A: Downloading Public Data (SRA)

If you are pulling public accessions down directly from the NCBI Sequence Read Archive (SRA):

```bash
# Download target SRA files
prefetch -O ./raw/ SRRxxxxxxxxx

# Convert SRA format into standard forward/reverse FASTQ files
fasterq-dump --split-files --threads 12 ./raw/SRRxxxxxxxxx -O ./fastq/

# Compress immediately using pigz (Parallel GZIP) to save ~70% disk space
pigz -p 12 ./fastq/SRRxxxxxxxxx_1.fastq ./fastq/SRRxxxxxxxxx_2.fastq
rm -rf ./raw/SRRxxxxxxxxx

```

### Mode B: Scaling Across Multiple Samples (Batch Scripting)

If you have a file with a list of multiple target IDs, do not run them manually one by one. Use a loop structure:

```bash
# Multi-sample execution loop example
for SAMPLE in $(cat samples.txt); do
    echo "Processing download & compression for: ${SAMPLE}"
    prefetch -O ./raw/ ${SAMPLE}
    fasterq-dump --split-files --threads 12 ./raw/${SAMPLE} -O ./fastq/
    pigz -p 12 ./fastq/${SAMPLE}_1.fastq ./fastq/${SAMPLE}_2.fastq
    rm -rf ./raw/${SAMPLE}
done

```

## 5. Quality Control & Adapter Trimming

Raw sequencing files contain technical artifacts like sequencing adapters or low-quality base calls near the end of the reads. We use `fastp` because it acts as an all-in-one quality filter and adapter trimmer.

```bash
fastp -i ./fastq/SRRxxxxxxxxx_1.fastq.gz -I ./fastq/SRRxxxxxxxxx_2.fastq.gz \
  -o ./qc/SRRxxxxxxxxx_R1.fastq.gz -O ./qc/SRRxxxxxxxxx_R2.fastq.gz \
  -h ./qc/SRRxxxxxxxxx_fastp.html -j ./qc/SRRxxxxxxxxx_fastp.json \
  --detect_adapter_for_pe --thread 12 --compression 6

```

> **What to check:** Open the generated `.html` file in your browser. Look for a clean drop in adapter content and ensure your mean quality score (Q30) stays high across the read length.

## 6. Alignment to the Reference Genome

This step matches your sequence fragments to their correct physical location on the human chromosomes.

> Standard pipelines output a massive, text-based intermediate `.sam` file, write it to the hard drive, and then read it back to convert it to a compressed binary `.bam` file. This wastes time and disk space. Instead, we use the Unix pipe operator (`|`) to stream the alignment output directly through the compression and sorting engines without ever writing a `.sam` file to disk.

```bash
bowtie2 -p 12 --very-sensitive -X 2000 -x ./refs/GRCh38_index \
  -1 ./qc/SRRxxxxxxxxx_R1.fastq.gz -2 ./qc/SRRxxxxxxxxx_R2.fastq.gz | \
  samtools view -@ 12 -bS - | \
  samtools sort -@ 12 -m 1G -T ./tmp/ -o ./bam/SRRxxxxxxxxx.aln.bam -

samtools index -@ 12 ./bam/SRRxxxxxxxxx.aln.bam

```

### Flag Breakdown:

* `--very-sensitive`: Maximizes alignment accuracy, ensuring we do not misassign reads in open regions.
* `-X 2000`: **Mandatory for ATAC-seq.** The default maximum insert size for Bowtie2 is 500bp. Because ATAC-seq cuts between nucleosomes, fragments wrapping around multiple nucleosomes can easily exceed 500bp. Setting this to 2000 ensures you do not artificially truncate your multi-nucleosome data.

## 7. Filtration, Deduplication, and the Tn5 Shift

ATAC-seq libraries suffer from specific background noise: mitochondrial DNA contamination (which lacks protective chromatin and is highly accessible to the Tn5 enzyme), low-quality alignments, PCR duplication artifacts, and technical blacklist zones.

```bash
# Step 1: Filter low-quality/unmapped reads and strip blacklist regions in a single stream
samtools view -@ 12 -u -q 30 -F 1804 ./bam/SRRxxxxxxxxx.aln.bam | \
  bedtools intersect -v -ubam -a - -b ./refs/ENCFF356LFX.bed | \

# Step 2: Collate, assign mate tags, and coordinate-sort for duplicate marking
  samtools collate -O -u -@ 12 - ./tmp/collate_prefix | \
  samtools fixmate -m -u -@ 12 - - | \
  samtools sort -@ 12 -T ./tmp/ -o ./bam/SRRxxxxxxxxx.fixmate.bam -

# Step 3: Remove PCR duplicates generated during library construction
samtools markdup -@ 12 -r ./bam/SRRxxxxxxxxx.fixmate.bam ./bam/SRRxxxxxxxxx.rmdup.bam
samtools index -@ 12 ./bam/SRRxxxxxxxxx.rmdup.bam

# Step 4: Apply the Tn5 Shift and re-sort prior to indexing (Strictly for downstream visualization tracks)
alignmentSieve -b ./bam/SRRxxxxxxxxx.rmdup.bam -o ./tmp/SRRxxxxxxxxx.shifted.unsorted.bam --ATACshift --numberOfProcessors 12
samtools sort -@ 12 -T ./tmp/ -o ./bam/SRRxxxxxxxxx.shifted.bam ./tmp/SRRxxxxxxxxx.shifted.unsorted.bam
samtools index -@ 12 ./bam/SRRxxxxxxxxx.shifted.bam

```


### Why do we shift?

The physical Tn5 transposase binds to open DNA as a dimer and inserts adapters with a characteristic **9-bp offset** (+4 bp on the forward strand, -5 bp on the reverse strand).

To see exactly where the enzyme cut at base-pair resolution on a genome browser, we use `alignmentSieve` to calculate a **shifted BAM**. However, **do not call peaks on this shifted BAM**. Modern peak callers like MACS3 are designed to read unshifted, paired-end configurations (`rmdup.bam`) and calculate the biological center point mathematically on their own.

## 8. Discovering Open Regions: Peak Calling

Now we identify areas where reads significantly cluster together compared to the genomic background. These clusters represent regions of open, accessible chromatin.

```bash
mkdir -p ./peaks/SRRxxxxxxxxx

macs3 callpeak -t ./bam/SRRxxxxxxxxx.rmdup.bam -f BAMPE -g hs \
  -n SRRxxxxxxxxx --outdir ./peaks/SRRxxxxxxxxx \
  --nomodel --shift -75 --extsize 150 --keep-dup all --qvalue 0.05

```

### Critical Parameter Insight:

* `-f BAMPE`: Instructs MACS3 to read the exact, real insertion lengths from your paired-end data rather than trying to guess the fragment sizes using single-end linear models.
* `--nomodel --shift -75 --extsize 150`: Bypasses the default internal shift models (which were originally designed for ChIP-seq data) and expands the signal around the cut site to properly resolve open chromatin dynamics.

## 9. Quantitative Signal Tracks (BigWigs)

To view your data in browsers like IGV or UCSC, you must convert your BAM alignments into continuous signal tracks. We use CPM (Counts Per Million mapped reads) normalization so that samples with different sequencing depths can be compared visually.

```bash
# Generate normalized continuous coverage tracks utilizing multi-threading
bamCoverage -b ./bam/SRRxxxxxxxxx.shifted.bam -o ./bw/SRRxxxxxxxxx.cpm.bw \
  --binSize 10 --normalizeUsing CPM --effectiveGenomeSize 2913022398 -p 12

```

* **`--binSize 10`**: Computes sequencing depth across consecutive 10-base-pair genomic bins. This prevents computational bloating while preserving sharp, nucleotide-level resolution of open promoter or enhancer hubs.
* **`--normalizeUsing CPM`**: Scales coverage scores to Counts Per Million mapped reads ($CPM = \frac{\text{reads in bin} \times 1,000,000}{\text{total mapped reads}}$). Without this normalization, a sample sequenced at 50 million reads would artificially appear to have twice as much chromatin openness as a sample sequenced at 25 million reads, even if their biology was identical.

## 10. Quality Assurance: Calculating your FRiP Score

How do you know if your experiment worked or if you just sequenced random genomic noise? The **FRiP score** (Fraction of Reads in Peaks) tells you what percentage of your total sequenced fragments actually landed inside called peaks.

```bash
# 1. Total clean, high-quality uniquely mapped reads (The Denominator)
TOTAL=$(samtools view -c -F 1804 -q 30 ./bam/SRRxxxxxxxxx.shifted.bam)

# 2. Total reads successfully intersecting with peak zones (The Numerator)
IN_PEAKS=$(bedtools intersect -u -abam ./bam/SRRxxxxxxxxx.shifted.bam -b ./peaks/SRRxxxxxxxxx/SRRxxxxxxxxx_peaks.narrowPeak | wc -l)

# 3. Calculate and print the percentage proportion
echo "Sample Quality FRiP Score:" $(echo "scale=4; $IN_PEAKS / $TOTAL" | bc)

```

| FRiP Range | Quality Assessment | Next Steps |
| --- | --- | --- |
| **> 0.20 (20%)** | **Excellent** | Exceptional library; proceed with confidence to differential analysis. |
| **0.15 - 0.20** | **Acceptable** | Standard quality; minor background noise present. |
| **< 0.15** | **Poor / Failed** | High background noise. Check for low cell viability or over-fragmentation during sample prep. |

## 11. Intermediate Cleanup

Once your final shifted BAM tracks, BigWigs, and peak lists are verified, clean up your workspace to free up disk space.

```bash
rm ./bam/*.aln.bam ./bam/*.filt.bam ./bam/*.clean.bam
rm -rf ./tmp/*

```

## Common Tradeoffs & Pitfalls to Remember

* **Replicates Matter:** While you can test a pipeline on a couple of samples, you cannot draw solid biological conclusions without a sufficient number of replicates. Always use at least three biological replicates per experimental group to ensure statistical confidence during downstream differential accessibility testing (e.g., using `DESeq2` or `DiffBind`).
* **Mitochondrial Reads:** Depending on the tissue type and preparation protocol, mitochondrial reads can sometimes make up over 50% of your raw data. Keep an eye on your read counts before and after the filtering step (`samtools view -F 1804`) to see how much of your library was lost to mitochondrial cleanup.
