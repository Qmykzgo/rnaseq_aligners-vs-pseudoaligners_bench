#!/bin/bash
# ==============================================================================
# 01_prepare_data.sh
# Purpose : Download GRCh38 (Ensembl 110) + ERCC spike-ins, extract the
#           transcriptome FASTA, build a chr21-masked GTF for Aim 3, and
#           download all benchmark FASTQ datasets from NCBI SRA.
#
# Reference:  GRCh38 primary assembly + Ensembl 110 GTF + ERCC-92 spike-ins
#
# Datasets:
#   Aim 1 — SEQC/MAQC-III  (GSE47774 / SRP025982)
#           SRR896663   Sample A (Universal HRR), BGI Rep1
#           SRR896679   Sample A (Universal HRR), BGI Rep2
#           SRR896743   Sample B (Brain Ref RNA),  BGI Rep1
#           SRR896759   Sample B (Brain Ref RNA),  BGI Rep2
#
#   Aim 2A — RNAontheBENCH iPSC lines (GSE63055)
#           SRR1643426  WBS phenotype Rep1
#           SRR1643427  WBS phenotype Rep2
#           SRR1643434  7dup phenotype Rep1
#
#   Aim 2B — MAQC qPCR validation (SRA010153)
#           Uses Aim 1 samples (SRR896663 / SRR896743) for downstream metrics
#
#   Aim 3  — Veeneman et al. 2015 deeply sequenced transcriptome
#           SRR1293902  Deeply sequenced polyA+ RNA (>60 M read pairs)
#
# Usage: bash 01_prepare_data.sh
# ==============================================================================
set -euo pipefail

CONDA_BASE="$(conda info --base 2>/dev/null || echo "${HOME}/miniconda3")"
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate rnaseq_bench_env

THREADS=8

# ---------- Directory layout ----------
BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="${BASE_DIR}/data"
REF_DIR="${DATA_DIR}/reference"
RAW_DIR="${DATA_DIR}/raw_fastqs"
SRA_CACHE="${DATA_DIR}/sra_cache"
TRUTH_DIR="${DATA_DIR}/truth_tables"

mkdir -p "${REF_DIR}" "${RAW_DIR}" "${SRA_CACHE}" "${TRUTH_DIR}"

# ==============================================================================
# 1.  REFERENCE GENOME  (GRCh38 Ensembl 110)
# ==============================================================================
echo "=== [1/6] Downloading GRCh38 Ensembl 110 FASTA + GTF ==="

GENOME_RAW="${REF_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
GTF_RAW="${REF_DIR}/Homo_sapiens.GRCh38.110.gtf"

if [ ! -f "${GENOME_RAW}" ]; then
    wget -c --progress=dot:giga \
        -O "${GENOME_RAW}.gz" \
        "http://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
    gunzip "${GENOME_RAW}.gz"
    echo "[OK] Genome FASTA decompressed."
else
    echo "[SKIP] Genome FASTA already present."
fi

if [ ! -f "${GTF_RAW}" ]; then
    wget -c --progress=dot:giga \
        -O "${GTF_RAW}.gz" \
        "http://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz"
    gunzip "${GTF_RAW}.gz"
    echo "[OK] GTF decompressed."
else
    echo "[SKIP] GTF already present."
fi

# ==============================================================================
# 2.  ERCC SPIKE-IN SEQUENCES (ERCC-92, ThermoFisher)
# ==============================================================================
echo "=== [2/6] Downloading and Integrating ERCC-92 Spike-ins ==="

ERCC_DIR="${REF_DIR}/ERCC92"
ERCC_URL="https://assets.thermofisher.com/TFS-Assets/LSG/manuals/ERCC92.zip"

if [ ! -f "${ERCC_DIR}/ERCC92.fa" ]; then
    wget -c --progress=dot:mega -O "${REF_DIR}/ERCC92.zip" "${ERCC_URL}"
    unzip -o "${REF_DIR}/ERCC92.zip" -d "${ERCC_DIR}"
    echo "[OK] ERCC92 extracted to ${ERCC_DIR}."
else
    echo "[SKIP] ERCC92 sequences already present."
fi

# Concatenate ERCC sequences into genome and annotation
GENOME_ERCC="${REF_DIR}/GRCh38_ERCC92.fa"
GTF_ERCC="${REF_DIR}/GRCh38_ERCC92.gtf"

if [ ! -f "${GENOME_ERCC}" ]; then
    cat "${GENOME_RAW}" "${ERCC_DIR}/ERCC92.fa" > "${GENOME_ERCC}"
    echo "[OK] ERCC-appended genome : ${GENOME_ERCC}"
fi

if [ ! -f "${GTF_ERCC}" ]; then
    cat "${GTF_RAW}" "${ERCC_DIR}/ERCC92.gtf" > "${GTF_ERCC}"
    echo "[OK] ERCC-appended GTF   : ${GTF_ERCC}"
fi

# ==============================================================================
# 3.  TRANSCRIPTOME FASTA  (for Salmon + kallisto)
# ==============================================================================
echo "=== [3/6] Extracting Transcriptome FASTA ==="

TRANSCRIPTOME="${REF_DIR}/GRCh38_ERCC92_transcriptome.fa"
if [ ! -f "${TRANSCRIPTOME}" ]; then
    # gffread writes spliced transcript sequences using the genome as template
    gffread "${GTF_ERCC}" -g "${GENOME_ERCC}" -w "${TRANSCRIPTOME}"
    echo "[OK] Transcriptome FASTA: ${TRANSCRIPTOME}"
else
    echo "[SKIP] Transcriptome FASTA already exists."
fi

# Build transcript-to-gene mapping (used by 04_aggregate_metrics.py)
TX2GENE="${REF_DIR}/tx2gene.tsv"
if [ ! -f "${TX2GENE}" ]; then
    awk 'BEGIN{OFS="\t"}
         $3 == "transcript" {
             tx=""; gn="";
             for(i=9;i<=NF;i++){
                 if($i=="transcript_id") tx=$(i+1);
                 if($i=="gene_id")       gn=$(i+1);
             }
             gsub(/[";]/,"",tx); gsub(/[";]/,"",gn);
             if(tx!="" && gn!="") print tx, gn
         }' "${GTF_ERCC}" > "${TX2GENE}"
    # ERCC entries: transcript == gene
    grep "^>" "${ERCC_DIR}/ERCC92.fa" | sed 's/^>//' \
        | awk '{split($0,a," "); print a[1]"\t"a[1]}' >> "${TX2GENE}"
    echo "[OK] tx2gene mapping written: ${TX2GENE}"
fi

# ==============================================================================
# 4.  CHR21-MASKED GTF  (Aim 3 novel-junction benchmark)
# ==============================================================================
# Strategy: Build alignment indices with annotations that EXCLUDE chromosome 21.
# Any reads originating from chr21 splice junctions will therefore appear as
# "novel" junctions, giving us a recall benchmark against a known truth set.
# In Ensembl GTF, chromosomes use bare numerals ("21", not "chr21").
echo "=== [4/6] Creating Chr21-Masked GTF for Aim 3 ==="

GTF_NO_CHR21="${REF_DIR}/GRCh38_ERCC92_noChr21.gtf"
if [ ! -f "${GTF_NO_CHR21}" ]; then
    # Keep comment lines (#) and all features NOT on chromosome 21
    awk '$0 ~ /^#/ || $1 != "21"' "${GTF_ERCC}" > "${GTF_NO_CHR21}"
    REMAINING=$(grep -c "^21	" "${GTF_NO_CHR21}" || true)
    echo "[OK] Chr21-masked GTF: ${GTF_NO_CHR21}  (chr21 lines remaining: ${REMAINING})"
fi

# Build chr21 truth junction set from the full GTF (used for precision/recall)
CHR21_JUNCTIONS="${REF_DIR}/chr21_true_junctions.tsv"
if [ ! -f "${CHR21_JUNCTIONS}" ]; then
    # Extract all exon pairs per transcript on chr21 → infer intron coordinates
    awk 'BEGIN{OFS="\t"}
         $1=="21" && $3=="exon" {
             tx="";
             for(i=9;i<=NF;i++){
                 if($i=="transcript_id"){tx=$(i+1); gsub(/[";]/,"",tx); break}
             }
             print tx, $1, $4, $5, $7
         }' "${GTF_ERCC}" \
    | sort -k1,1 -k3,3n \
    | awk 'BEGIN{OFS="\t"; prev_tx=""; prev_end=0; prev_chr=""; prev_strand=""}
           {
               tx=$1; chr=$2; start=$3; end=$4; strand=$5;
               if(tx == prev_tx && chr == prev_chr){
                   # intron runs from prev_end+1 to start-1
                   print chr, prev_end+1, start-1, strand
               }
               prev_tx=tx; prev_end=end; prev_chr=chr; prev_strand=strand
           }' \
    | sort -k1,1 -k2,2n -u \
    > "${CHR21_JUNCTIONS}"
    echo "[OK] Chr21 truth junctions: ${CHR21_JUNCTIONS}  ($(wc -l < "${CHR21_JUNCTIONS}") entries)"
fi

# ==============================================================================
# 5.  DOWNLOAD BENCHMARK FASTQs FROM NCBI SRA
# ==============================================================================
echo "=== [5/6] Downloading SRA Datasets ==="

# All accessions: associative arrays (accession → human-readable label)
declare -A AIM1=( [SRR896663]="SampleA_BGI_Rep1"
                  [SRR896679]="SampleA_BGI_Rep2"
                  [SRR896743]="SampleB_BGI_Rep1"
                  [SRR896759]="SampleB_BGI_Rep2" )

declare -A AIM2A=( [SRR1643426]="WBS_Rep1"
                   [SRR1643427]="WBS_Rep2"
                   [SRR1643434]="7dup_Rep1" )

declare -A AIM2B=( )

declare -A AIM3=( [SRR1293902]="Veeneman_deepseq" )

# Generic download function
download_sra() {
    local SRR="$1"
    local LABEL="$2"

    # Check if already done (both paired-end files expected)
    if [ -f "${RAW_DIR}/${SRR}_1.fastq.gz" ]; then
        echo "[SKIP] ${SRR} (${LABEL}) — already downloaded."
        return 0
    fi

    if [ "${SRR}" = "SRR1293902" ]; then
        echo "[INFO] Downloading subset (100k spots) of ${SRR} (${LABEL}) directly ..."
        fastq-dump -X 100000 \
            -O "${RAW_DIR}" \
            --split-files \
            --gzip \
            "${SRR}"
        echo "[OK] ${SRR} (${LABEL}) subset — done."
        return 0
    fi

    echo "[INFO] Prefetching ${SRR} (${LABEL}) ..."
    prefetch "${SRR}" \
        -O "${SRA_CACHE}" \
        --max-size 60G \
        --progress

    echo "[INFO] Converting ${SRR} → FASTQ ..."
    fasterq-dump "${SRA_CACHE}/${SRR}" \
        -O "${RAW_DIR}" \
        --threads "${THREADS}" \
        --split-files \
        --progress

    echo "[INFO] Compressing ${SRR} FASTQs ..."
    # Handle both paired-end and single-end outputs robustly
    for f in "${RAW_DIR}/${SRR}"*.fastq; do
        [ -f "${f}" ] && gzip -f "${f}"
    done

    echo "[OK] ${SRR} (${LABEL}) — done."
}

echo "--- Aim 1: SEQC/MAQC-III (4 samples) ---"
for SRR in "${!AIM1[@]}"; do download_sra "${SRR}" "${AIM1[$SRR]}"; done

echo "--- Aim 2A: RNAontheBENCH (3 samples) ---"
for SRR in "${!AIM2A[@]}"; do download_sra "${SRR}" "${AIM2A[$SRR]}"; done

echo "--- Aim 2B: MAQC qPCR validation (2 samples) ---"
for SRR in "${!AIM2B[@]}"; do download_sra "${SRR}" "${AIM2B[$SRR]}"; done

echo "--- Aim 3: Veeneman deep-sequencing (1 sample) ---"
for SRR in "${!AIM3[@]}"; do download_sra "${SRR}" "${AIM3[$SRR]}"; done

# ==============================================================================
# 6.  DOWNLOAD ACCURACY TRUTH TABLES
# ==============================================================================
echo "=== [6/6] Downloading Accuracy Truth Tables ==="

# Aim 2B: Everaert et al. 2017 whole-transcriptome qPCR
# GEO accession: GSE83402 — normalized Cq values for 18,080 protein-coding genes.
QPCR_FILE="${TRUTH_DIR}/GSE83402_qPCR_normalized.txt"
if [ ! -f "${QPCR_FILE}" ]; then
    QPCR_GEO_URL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE83nnn/GSE83402/suppl/GSE83402_MAQCA_MAQCB_qPCR_RefSeq_normalized.txt.gz"
    wget -c -O "${QPCR_FILE}.gz" "${QPCR_GEO_URL}" \
        && gunzip "${QPCR_FILE}.gz" \
        && echo "[OK] qPCR truth table downloaded." \
        || echo "[WARN] Auto-download failed. Manually download from:
           https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE83402
           Save the normalized Cq table as: ${QPCR_FILE}"
else
    echo "[SKIP] qPCR truth table already present."
fi

# Aim 2A: RNAontheBENCH NanoString truth (CSV from GitHub)
BENCH_FILE="${TRUTH_DIR}/RNAontheBENCH_truth_counts.csv"
if [ ! -f "${BENCH_FILE}" ]; then
    BENCH_URL="https://raw.githubusercontent.com/plger/RNAontheBENCH/master/data/truth.counts.csv"
    wget -c -O "${BENCH_FILE}" "${BENCH_URL}" \
        && echo "[OK] RNAontheBENCH truth counts downloaded." \
        || echo "[WARN] Auto-download failed. Clone https://github.com/plger/RNAontheBENCH
           and copy data/truth.counts.csv to: ${BENCH_FILE}"
else
    echo "[SKIP] RNAontheBENCH truth already present."
fi

# ==============================================================================
# SUMMARY
# ==============================================================================
echo ""
echo "============================================================"
echo "  [DONE] Data preparation complete."
echo "============================================================"
printf "  %-30s %s\n" "Genome + ERCC :"  "${GENOME_ERCC}"
printf "  %-30s %s\n" "GTF + ERCC :"     "${GTF_ERCC}"
printf "  %-30s %s\n" "GTF no-chr21 :"   "${GTF_NO_CHR21}"
printf "  %-30s %s\n" "Chr21 junctions:" "${CHR21_JUNCTIONS}"
printf "  %-30s %s\n" "Transcriptome :"  "${TRANSCRIPTOME}"
printf "  %-30s %s\n" "tx2gene map :"    "${TX2GENE}"
printf "  %-30s %s\n" "FASTQs :"         "${RAW_DIR}"
printf "  %-30s %s\n" "Truth tables :"   "${TRUTH_DIR}"
