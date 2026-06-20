#!/bin/bash
# ==============================================================================
# 02_run_alignment.sh
# Purpose : Build STAR and HISAT2 indices, run alignments for all datasets,
#           and quantify with featureCounts.
#
# Index strategy:
#   star_idx/          — Full ERCC-appended genome, used for Aim 1 & 2
#   star_idx_noChr21/  — Chr21-masked annotation, used for Aim 3 (novel junctions)
#   hisat2_idx/        — Full ERCC-appended genome, used for Aim 1 & 2
#   hisat2_idx_noChr21/— Chr21-masked splice sites, used for Aim 3
#
# All compute-heavy commands are wrapped with /usr/bin/time -v to capture
# wall-clock time and peak RAM to dedicated .log files under results/logs/.
#
# Usage: bash 02_run_alignment.sh
# ==============================================================================
set -euo pipefail

CONDA_BASE="$(conda info --base 2>/dev/null || echo "${HOME}/miniconda3")"
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate rnaseq_bench_env

THREADS=8

BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="${BASE_DIR}/data"
REF_DIR="${DATA_DIR}/reference"
RAW_DIR="${DATA_DIR}/raw_fastqs"
OUT_DIR="${BASE_DIR}/results/alignment"
LOG_DIR="${BASE_DIR}/results/logs"

mkdir -p "${OUT_DIR}/star" "${OUT_DIR}/hisat2" "${LOG_DIR}"

echo "=== Decompressing FASTQs if necessary ==="
for f in "${RAW_DIR}"/*.fastq.gz; do
    if [ -f "$f" ]; then
        echo "Decompressing $f ..."
        gunzip -f "$f"
    fi
done

GENOME_ERCC="${REF_DIR}/GRCh38_ERCC92.fa"
GTF_ERCC="${REF_DIR}/GRCh38_ERCC92.gtf"
GTF_NO_CHR21="${REF_DIR}/GRCh38_ERCC92_noChr21.gtf"

# sjdbOverhang = max_read_length - 1 (use 99 for 100 bp SEQC reads; 149 for 150 bp)
SJDB_OVERHANG=99

# ==============================================================================
# HELPER: detect read length from first FASTQ
# ==============================================================================
get_read_length() {
    local FQ="$1"
    # Read the second line (sequence) of the first read; strip whitespace
    if [[ "${FQ}" == *.gz ]]; then
        set +o pipefail
        local LEN
        LEN=$(zcat "${FQ}" | awk 'NR==2{print length($0); exit}')
        set -o pipefail
        echo "${LEN}"
    else
        awk 'NR==2{print length($0); exit}' "${FQ}"
    fi
}

# ==============================================================================
# 1.  INDEX BUILDING
# ==============================================================================

# ---------- 1a. HISAT2 full index (Aim 1 & 2) ----------
HISAT2_IDX="${REF_DIR}/hisat2_idx/genome"
if [ ! -f "${HISAT2_IDX}.1.ht2" ]; then
    echo "=== Building HISAT2 full index ==="
    mkdir -p "${REF_DIR}/hisat2_idx"

    # Extract known splice sites and exon boundaries for the genome-directed index
    hisat2_extract_splice_sites.py "${GTF_ERCC}" \
        > "${REF_DIR}/hisat2_idx/splice_sites.txt"
    hisat2_extract_exons.py "${GTF_ERCC}" \
        > "${REF_DIR}/hisat2_idx/exons.txt"

    /usr/bin/time -v \
        -o "${LOG_DIR}/hisat2_index_full.log" \
        hisat2-build \
            -p "${THREADS}" \
            "${GENOME_ERCC}" \
            "${HISAT2_IDX}" \
        2>&1 | tee "${LOG_DIR}/hisat2_index_full.stdout.log"
    echo "[OK] HISAT2 full index done."
else
    echo "[SKIP] HISAT2 full index already built."
fi

# ---------- 1b. HISAT2 chr21-masked index (Aim 3) ----------
HISAT2_IDX_NC21="${REF_DIR}/hisat2_idx_noChr21/genome"
if [ ! -f "${HISAT2_IDX_NC21}.1.ht2" ]; then
    echo "=== Building HISAT2 chr21-masked index (Aim 3) ==="
    mkdir -p "${REF_DIR}/hisat2_idx_noChr21"

    hisat2_extract_splice_sites.py "${GTF_NO_CHR21}" \
        > "${REF_DIR}/hisat2_idx_noChr21/splice_sites.txt"
    hisat2_extract_exons.py "${GTF_NO_CHR21}" \
        > "${REF_DIR}/hisat2_idx_noChr21/exons.txt"

    /usr/bin/time -v \
        -o "${LOG_DIR}/hisat2_index_noChr21.log" \
        hisat2-build \
            -p "${THREADS}" \
            "${GENOME_ERCC}" \
            "${HISAT2_IDX_NC21}" \
        2>&1 | tee "${LOG_DIR}/hisat2_index_noChr21.stdout.log"
    echo "[OK] HISAT2 chr21-masked index done."
else
    echo "[SKIP] HISAT2 chr21-masked index already built."
fi

# ---------- 1c. STAR full index (Aim 1 & 2) ----------
STAR_IDX="${REF_DIR}/star_idx"
if [ ! -f "${STAR_IDX}/SA" ]; then
    echo "=== Building STAR full index ==="
    mkdir -p "${STAR_IDX}"
    /usr/bin/time -v \
        -o "${LOG_DIR}/star_index_full.log" \
        STAR \
            --runThreadN "${THREADS}" \
            --runMode genomeGenerate \
            --genomeDir "${STAR_IDX}" \
            --genomeFastaFiles "${GENOME_ERCC}" \
            --sjdbGTFfile "${GTF_ERCC}" \
            --sjdbOverhang "${SJDB_OVERHANG}" \
            --limitGenomeGenerateRAM 10000000000 \
            --genomeChrBinNbits 16 \
            --genomeSAindexNbases 11 \
        2>&1 | tee "${LOG_DIR}/star_index_full.stdout.log"
    echo "[OK] STAR full index done."
else
    echo "[SKIP] STAR full index already built."
fi

# ---------- 1d. STAR chr21-masked index (Aim 3) ----------
STAR_IDX_NC21="${REF_DIR}/star_idx_noChr21"
if [ ! -f "${STAR_IDX_NC21}/SA" ]; then
    echo "=== Building STAR chr21-masked index (Aim 3) ==="
    mkdir -p "${STAR_IDX_NC21}"
    /usr/bin/time -v \
        -o "${LOG_DIR}/star_index_noChr21.log" \
        STAR \
            --runThreadN "${THREADS}" \
            --runMode genomeGenerate \
            --genomeDir "${STAR_IDX_NC21}" \
            --genomeFastaFiles "${GENOME_ERCC}" \
            --sjdbGTFfile "${GTF_NO_CHR21}" \
            --sjdbOverhang "${SJDB_OVERHANG}" \
            --limitGenomeGenerateRAM 10000000000 \
            --genomeChrBinNbits 16 \
            --genomeSAindexNbases 11 \
        2>&1 | tee "${LOG_DIR}/star_index_noChr21.stdout.log"
    echo "[OK] STAR chr21-masked index done."
else
    echo "[SKIP] STAR chr21-masked index already built."
fi

# ==============================================================================
# 2.  ALIGNMENT FUNCTION  (Aim 1 & 2 — standard one-pass alignment)
# ==============================================================================
run_alignment_standard() {
    local SAMPLE="$1"   # SRR accession
    local R1="${RAW_DIR}/${SAMPLE}_1.fastq"
    local R2="${RAW_DIR}/${SAMPLE}_2.fastq"

    # Require both FASTQ files
    if [ ! -f "${R1}" ] || [ ! -f "${R2}" ]; then
        echo "[WARN] FASTQs not found for ${SAMPLE} — skipping."
        return 1
    fi

    mkdir -p "${OUT_DIR}/star/${SAMPLE}" "${OUT_DIR}/hisat2"

    # ---- STAR 1-pass ----
    if [ ! -f "${OUT_DIR}/star/${SAMPLE}/Aligned.sortedByCoord.out.bam" ]; then
        echo "[RUN] STAR alignment — ${SAMPLE}"
        RL=$(get_read_length "${R1}")
        echo "      Detected read length = ${RL} bp  (sjdbOverhang set to $((RL-1)))"
        /usr/bin/time -v \
            -o "${LOG_DIR}/star_align_${SAMPLE}.log" \
            STAR \
                --runThreadN "${THREADS}" \
                --genomeDir "${STAR_IDX}" \
                --readFilesIn "${R1}" "${R2}" \
                --outSAMtype BAM Unsorted \
                --outSAMattributes NH HI AS NM MD \
                --outFileNamePrefix "${OUT_DIR}/star/${SAMPLE}/" \
                --quantMode GeneCounts \
            2>&1 | tee "${LOG_DIR}/star_align_${SAMPLE}.stdout.log"
        samtools sort -@ "${THREADS}" -o "${OUT_DIR}/star/${SAMPLE}/Aligned.sortedByCoord.out.bam" "${OUT_DIR}/star/${SAMPLE}/Aligned.out.bam"
        rm "${OUT_DIR}/star/${SAMPLE}/Aligned.out.bam"
        samtools index "${OUT_DIR}/star/${SAMPLE}/Aligned.sortedByCoord.out.bam"
        echo "[OK] STAR — ${SAMPLE}"
    else
        echo "[SKIP] STAR alignment already done for ${SAMPLE}."
    fi

    # ---- HISAT2 ----
    local HISAT2_BAM="${OUT_DIR}/hisat2/${SAMPLE}.bam"
    if [ ! -f "${HISAT2_BAM}" ]; then
        echo "[RUN] HISAT2 alignment — ${SAMPLE}"
        /usr/bin/time -v \
            -o "${LOG_DIR}/hisat2_align_${SAMPLE}.log" \
            hisat2 \
                -p "${THREADS}" \
                -x "${HISAT2_IDX}" \
                -1 "${R1}" -2 "${R2}" \
                --dta \
                --novel-splicesite-outfile \
                    "${OUT_DIR}/hisat2/${SAMPLE}_novel_junctions.txt" \
            2>"${LOG_DIR}/hisat2_align_${SAMPLE}.stdout.log" \
            | samtools sort \
                -@ "${THREADS}" \
                -o "${HISAT2_BAM}" -
        samtools index "${HISAT2_BAM}"
        echo "[OK] HISAT2 — ${SAMPLE}"
    else
        echo "[SKIP] HISAT2 alignment already done for ${SAMPLE}."
    fi

    # ---- featureCounts (on both BAMs simultaneously for this sample) ----
    local FC_STAR="${OUT_DIR}/star/${SAMPLE}_featureCounts.txt"
    if [ ! -f "${FC_STAR}" ]; then
        echo "[RUN] featureCounts (STAR BAM) — ${SAMPLE}"
        /usr/bin/time -v \
            -o "${LOG_DIR}/featureCounts_star_${SAMPLE}.log" \
            featureCounts \
                -T "${THREADS}" \
                -p \
                -B \
                -C \
                -s 0 \
                -a "${GTF_ERCC}" \
                -o "${FC_STAR}" \
                "${OUT_DIR}/star/${SAMPLE}/Aligned.sortedByCoord.out.bam" \
            2>&1 | tee "${LOG_DIR}/featureCounts_star_${SAMPLE}.stdout.log"
        echo "[OK] featureCounts STAR — ${SAMPLE}"
    else
        echo "[SKIP] featureCounts (STAR) already done for ${SAMPLE}."
    fi

    local FC_HISAT2="${OUT_DIR}/hisat2/${SAMPLE}_featureCounts.txt"
    if [ ! -f "${FC_HISAT2}" ]; then
        echo "[RUN] featureCounts (HISAT2 BAM) — ${SAMPLE}"
        /usr/bin/time -v \
            -o "${LOG_DIR}/featureCounts_hisat2_${SAMPLE}.log" \
            featureCounts \
                -T "${THREADS}" \
                -p \
                -B \
                -C \
                -s 0 \
                -a "${GTF_ERCC}" \
                -o "${FC_HISAT2}" \
                "${HISAT2_BAM}" \
            2>&1 | tee "${LOG_DIR}/featureCounts_hisat2_${SAMPLE}.stdout.log"
        echo "[OK] featureCounts HISAT2 — ${SAMPLE}"
    else
        echo "[SKIP] featureCounts (HISAT2) already done for ${SAMPLE}."
    fi
}

# ==============================================================================
# 3.  RUN STANDARD ALIGNMENT  (Aim 1 & 2 samples)
# ==============================================================================
AIM1_SAMPLES=( SRR896663 SRR896679 SRR896743 SRR896759 )
AIM2A_SAMPLES=( SRR1643426 SRR1643427 SRR1643434 )
AIM2B_SAMPLES=( )

ALL_STD_SAMPLES=( "${AIM1_SAMPLES[@]}" "${AIM2A_SAMPLES[@]}" "${AIM2B_SAMPLES[@]}" )

echo "=== Running standard alignment (Aim 1 & 2) ==="
for SAMPLE in "${ALL_STD_SAMPLES[@]}"; do
    run_alignment_standard "${SAMPLE}"
done

# ==============================================================================
# 4.  AIM 3 — STAR TWO-PASS + HISAT2  (chr21-masked index, junction benchmark)
# ==============================================================================
AIM3_SAMPLE="SRR1293902"
if [ ! -f "${RAW_DIR}/${AIM3_SAMPLE}_1.fastq" ]; then
    gunzip -c "${RAW_DIR}/${AIM3_SAMPLE}_1.fastq.gz" > "${RAW_DIR}/${AIM3_SAMPLE}_1.fastq"
    gunzip -c "${RAW_DIR}/${AIM3_SAMPLE}_2.fastq.gz" > "${RAW_DIR}/${AIM3_SAMPLE}_2.fastq"
fi

R1_A3="${RAW_DIR}/${AIM3_SAMPLE}_1.fastq"
R2_A3="${RAW_DIR}/${AIM3_SAMPLE}_2.fastq"

echo "=== Aim 3: STAR two-pass alignment (${AIM3_SAMPLE}) ==="

mkdir -p "${OUT_DIR}/star/${AIM3_SAMPLE}_2pass"

if [ ! -f "${OUT_DIR}/star/${AIM3_SAMPLE}_2pass/Aligned.sortedByCoord.out.bam" ]; then
    RL_A3=$(get_read_length "${R1_A3}")
    echo "  Detected read length = ${RL_A3} bp"
    /usr/bin/time -v \
        -o "${LOG_DIR}/star_2pass_align_${AIM3_SAMPLE}.log" \
        STAR \
            --runThreadN "${THREADS}" \
            --genomeDir "${STAR_IDX_NC21}" \
            --readFilesIn "${R1_A3}" "${R2_A3}" \
            --twopassMode Basic \
            --outSAMtype BAM Unsorted \
            --outSAMattributes NH HI AS NM MD \
            --outFileNamePrefix "${OUT_DIR}/star/${AIM3_SAMPLE}_2pass/" \
        2>&1 | tee "${LOG_DIR}/star_2pass_align_${AIM3_SAMPLE}.stdout.log"
    samtools sort -@ "${THREADS}" -o "${OUT_DIR}/star/${AIM3_SAMPLE}_2pass/Aligned.sortedByCoord.out.bam" "${OUT_DIR}/star/${AIM3_SAMPLE}_2pass/Aligned.out.bam"
    rm "${OUT_DIR}/star/${AIM3_SAMPLE}_2pass/Aligned.out.bam"
    samtools index "${OUT_DIR}/star/${AIM3_SAMPLE}_2pass/Aligned.sortedByCoord.out.bam"
    echo "[OK] STAR 2-pass — ${AIM3_SAMPLE}"
else
    echo "[SKIP] STAR 2-pass already done."
fi

# STAR 1-pass for the same sample (baseline comparison)
mkdir -p "${OUT_DIR}/star/${AIM3_SAMPLE}_1pass"
if [ ! -f "${OUT_DIR}/star/${AIM3_SAMPLE}_1pass/Aligned.sortedByCoord.out.bam" ]; then
    /usr/bin/time -v \
        -o "${LOG_DIR}/star_1pass_align_${AIM3_SAMPLE}.log" \
        STAR \
            --runThreadN "${THREADS}" \
            --genomeDir "${STAR_IDX_NC21}" \
            --readFilesIn "${R1_A3}" "${R2_A3}" \
            --outSAMtype BAM Unsorted \
            --outSAMattributes NH HI AS NM MD \
            --outFileNamePrefix "${OUT_DIR}/star/${AIM3_SAMPLE}_1pass/" \
    2>&1 | tee "${LOG_DIR}/star_1pass_align_${AIM3_SAMPLE}.stdout.log"
    samtools sort -@ "${THREADS}" -o "${OUT_DIR}/star/${AIM3_SAMPLE}_1pass/Aligned.sortedByCoord.out.bam" "${OUT_DIR}/star/${AIM3_SAMPLE}_1pass/Aligned.out.bam"
    rm "${OUT_DIR}/star/${AIM3_SAMPLE}_1pass/Aligned.out.bam"
    samtools index "${OUT_DIR}/star/${AIM3_SAMPLE}_1pass/Aligned.sortedByCoord.out.bam"
    echo "[OK] STAR 1-pass (baseline) — ${AIM3_SAMPLE}"
else
    echo "[SKIP] STAR 1-pass baseline already done."
fi

# HISAT2 with chr21-masked splice sites (Aim 3)
HISAT2_BAM_A3="${OUT_DIR}/hisat2/${AIM3_SAMPLE}_noChr21.bam"
if [ ! -f "${HISAT2_BAM_A3}" ]; then
    echo "[RUN] HISAT2 Aim 3 alignment — ${AIM3_SAMPLE}"
    /usr/bin/time -v \
        -o "${LOG_DIR}/hisat2_aim3_align_${AIM3_SAMPLE}.log" \
        hisat2 \
            -p "${THREADS}" \
            -x "${HISAT2_IDX_NC21}" \
            -1 "${R1_A3}" -2 "${R2_A3}" \
            --dta \
            --novel-splicesite-outfile \
                "${OUT_DIR}/hisat2/${AIM3_SAMPLE}_novel_junctions.txt" \
        2>"${LOG_DIR}/hisat2_aim3_align_${AIM3_SAMPLE}.stdout.log" \
        | samtools sort \
            -@ "${THREADS}" \
            -o "${HISAT2_BAM_A3}" -
    samtools index "${HISAT2_BAM_A3}"
    echo "[OK] HISAT2 Aim 3 — ${AIM3_SAMPLE}"
else
    echo "[SKIP] HISAT2 Aim 3 already done."
fi

echo ""
echo "[DONE] All alignments and featureCounts complete."
echo "  Results: ${OUT_DIR}"
echo "  Logs   : ${LOG_DIR}"
