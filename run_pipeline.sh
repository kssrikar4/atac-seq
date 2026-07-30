#!/bin/bash
set -e

printf "===============================================\n"
printf "       AUTOMATED BULK ATAC-SEQ PIPELINE        \n"
printf "===============================================\n\n"

if [ -z "$CONDA_PREFIX" ] || [[ "$CONDA_PREFIX" != */env ]]; then
    printf "[ERROR] Local Conda environment './env' is not active.\n"
    printf "Please run: conda activate ./env\n"
    exit 1
fi

printf "Enter SRA Accession Number(s) (separated by spaces if multiple): "
read -r -a SAMPLES

if [ ${#SAMPLES[@]} -eq 0 ]; then
    printf "[ERROR] No SRA accessions provided. Exiting.\n"
    exit 1
fi

printf "\n[1/11] Creating directory infrastructure...\n"
mkdir -p ./refs ./raw ./fastq ./qc ./bam ./peaks ./bw ./tmp

export TMPDIR="$(pwd)/tmp"
set -euo pipefail
ulimit -Sn 65536 2>/dev/null || ulimit -Sn "$(ulimit -Hn)"

if [ ! -f "./refs/GRCh38_index.1.bt2" ]; then
    printf "\n[2/11] Reference genome indices not found. Initiating downloads...\n"
    wget -P ./refs/ https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GCA_000001405.15_GRCh38_genomic.fna.gz
    gzip -d ./refs/GCA_000001405.15_GRCh38_genomic.fna.gz

    printf "Building Bowtie2 index (this may take up to 30 minutes)...\n"
    bowtie2-build --threads 12 ./refs/GCA_000001405.15_GRCh38_genomic.fna ./refs/GRCh38_index

    samtools faidx ./refs/GCA_000001405.15_GRCh38_genomic.fna
    cut -f1,2 ./refs/GCA_000001405.15_GRCh38_genomic.fna.fai > ./refs/GRCh38.chrom.sizes

    wget -P ./refs/ https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz
    gzip -d ./refs/ENCFF356LFX.bed.gz
else
    printf "\n[2/11] Reference structures and indices verified.\n"
fi

for SAMPLE in "${SAMPLES[@]}"; do
    printf "\n===============================================\n"
    printf "STARTING PROCESSING FOR SAMPLE: %s\n" "$SAMPLE"
    printf "===============================================\n"

    printf "\n[3/11] Fetching raw sequence files from SRA...\n"
    prefetch -O ./raw/ "$SAMPLE"
    fasterq-dump --split-files --threads 12 ./raw/"$SAMPLE" -O ./fastq/
    pigz -p 12 ./fastq/"${SAMPLE}"_1.fastq ./fastq/"${SAMPLE}"_2.fastq
    rm -rf ./raw/"$SAMPLE"

    printf "\n[4/11] Executing quality control and adapter trimming...\n"
    fastp -i ./fastq/"${SAMPLE}"_1.fastq.gz -I ./fastq/"${SAMPLE}"_2.fastq.gz \
      -o ./qc/"${SAMPLE}"_R1.fastq.gz -O ./qc/"${SAMPLE}"_R2.fastq.gz \
      -h ./qc/"${SAMPLE}"_fastp.html -j ./qc/"${SAMPLE}"_fastp.json \
      --detect_adapter_for_pe --thread 12 --compression 6

    printf "\n[5/11] Aligning reads to reference genome (piping directly to sort)...\n"
    bowtie2 -p 12 --very-sensitive -X 2000 -x ./refs/GRCh38_index \
      -1 ./qc/"${SAMPLE}"_R1.fastq.gz -2 ./qc/"${SAMPLE}"_R2.fastq.gz | \
      samtools view -@ 12 -bS - | \
      samtools sort -@ 12 -m 1G -T ./tmp/ -o ./bam/"${SAMPLE}".aln.bam -
    samtools index -@ 12 ./bam/"${SAMPLE}".aln.bam

    printf "\n[6/11] Filtering low quality reads and resolving ENCODE blacklist...\n"
    samtools view -@ 12 -b -q 30 -F 1804 ./bam/"${SAMPLE}".aln.bam > ./bam/"${SAMPLE}".filt.bam
    bedtools intersect -v -abam ./bam/"${SAMPLE}".filt.bam -b ./refs/ENCFF356LFX.bed | \
      samtools sort -@ 12 -T ./tmp/ -o ./bam/"${SAMPLE}".clean.bam -
    samtools index -@ 12 ./bam/"${SAMPLE}".clean.bam

    printf "\n[7/11] Fixing mate tags, marking PCR duplicates, and executing Tn5 shifting...\n"
    samtools collate -O -u -@ 12 ./bam/"${SAMPLE}".clean.bam ./tmp/collate_prefix_"${SAMPLE}" | \
      samtools fixmate -m -u -@ 12 - - | \
      samtools sort -@ 12 -T ./tmp/ -o ./bam/"${SAMPLE}".fixmate.bam -

    samtools markdup -@ 12 -r ./bam/"${SAMPLE}".fixmate.bam ./bam/"${SAMPLE}".rmdup.bam
    samtools index -@ 12 ./bam/"${SAMPLE}".rmdup.bam
    alignmentSieve -b ./bam/"${SAMPLE}".rmdup.bam -o ./bam/"${SAMPLE}".shifted.bam --ATACshift --numberOfProcessors 12
    samtools index -@ 12 ./bam/"${SAMPLE}".shifted.bam

    printf "\n[8/11] Calling open chromatin peaks via MACS3...\n"
    mkdir -p ./peaks/"$SAMPLE"
    macs3 callpeak -t ./bam/"${SAMPLE}".rmdup.bam -f BAMPE -g hs \
      -n "$SAMPLE" --outdir ./peaks/"$SAMPLE" \
      --nomodel --shift -75 --extsize 150 --keep-dup all --qvalue 0.05

    printf "\n[9/11] Generating normalized continuous signal tracks (BigWig)...\n"
    bamCoverage -b ./bam/"${SAMPLE}".shifted.bam -o ./bw/"${SAMPLE}".cpm.bw \
      --binSize 10 --normalizeUsing CPM --effectiveGenomeSize 2913022398 -p 12

    printf "\n[10/11] Calculating Fraction of Reads in Peaks (FRiP) metric...\n"
    TOTAL=$(samtools view -c -F 1804 -q 30 ./bam/"${SAMPLE}".shifted.bam)
    IN_PEAKS=$(bedtools intersect -u -abam ./bam/"${SAMPLE}".shifted.bam -b ./peaks/"${SAMPLE}"/"${SAMPLE}"_peaks.narrowPeak | wc -l)
    FRIP=$(echo "scale=4; $IN_PEAKS / $TOTAL" | bc)
    printf "%s FRiP Quality Metric: %s\n" "$SAMPLE" "$FRIP"
    echo "$FRIP" > ./qc/"${SAMPLE}"_frip.txt

done

printf "\n[11/11] Cleaning up intermediate binary alignments...\n"
rm -f ./bam/*.aln.bam ./bam/*.filt.bam ./bam/*.clean.bam
rm -rf ./tmp/*

printf "\nPipeline execution completed successfully for all samples.\n"
