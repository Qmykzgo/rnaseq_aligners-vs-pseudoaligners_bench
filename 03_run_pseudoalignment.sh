#!/bin/bash
# ==============================================================================
# 03_run_pseudoalignment.sh
# Purpose : Build Salmon (decoy-aware) and kallisto indices, then quantify
#           Aim 1 and Aim 2 samples with both tools.
#
# Note    : Pseudoaligners are not designed for novel-splice-junction discovery;
#           Aim 3 is intentionally excluded here and covered in 02_run_alignment.sh.
#
# Salmon decoy strategy: The full genome FASTA is appended to the transcriptome
#   and a decoys.txt list is provided. This prevents spurious alignment of reads
#   that map to non-transcribed genomic regions, improving TPM accuracy.
#
# All compute-heavy commands are wrapped with /usr/bin/time -v.
#
# Usage: bash 03_run_pseudoalignment.sh
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
OUT_DIR="${BASE_DIR}/results/pseudoalignment"
LOG_DIR="${BASE_DIR}/results/logs"

mkdir -p "${OUT_DIR}/salmon" "${OUT_DIR}/kallisto" "${LOG_DIR}"

GENOME_ERCC="${REF_DIR}/GRCh38_ERCC92.fa"
TRANSCRIPTOME="${REF_DIR}/GRCh38_ERCC92_transcriptome.fa"

SALMON_IDX="${REF_DIR}/salmon_idx"
KALLISTO_IDX="${REF_DIR}/kallisto_idx.idx"

BOOTSTRAP=100  # kallisto bootstrap replicates for uncertainty estimation

# ==============================================================================
# 1.  SALMON DECOY-AWARE INDEX
#     Step A: Build the "gentrome" — transcriptome + full genome concatenated.
#     Step B: Create decoys.txt listing all genome contig names.
#     Step C: Index with --decoys flag.
#     This is the recommended salmon indexing approach (Srivastava et al., 2020).
# ==============================================================================
if [ ! -d "${SALMON_IDX}" ]; then
    echo "=== Building Salmon decoy-aware index ==="

    GENTROME="${REF_DIR}/gentrome.fa"
    DECOYS="${REF_DIR}/decoys.txt"

    # Extract only chromosome 21 sequence (the decoy) from GENOME_ERCC
    awk '/^>ERCC/{exit} {print}' "${GENOME_ERCC}" > "${REF_DIR}/chr21_only.fa"
    echo "21" > "${DECOYS}"
    echo "[INFO] Decoy contigs: $(cat "${DECOYS}")"

    # Concatenate transcriptome then chromosome 21 sequence (decoy is contiguous at the end)
    cat "${TRANSCRIPTOME}" "${REF_DIR}/chr21_only.fa" > "${GENTROME}"
    echo "[INFO] Gentrome size: $(du -sh "${GENTROME}" | cut -f1)"

    /usr/bin/time -v \
        -o "${LOG_DIR}/salmon_index.log" \
        salmon index \
            -t "${GENTROME}" \
            -d "${DECOYS}" \
            -i "${SALMON_IDX}" \
            --threads "${THREADS}" \
            --gencode \
        2>&1 | tee "${LOG_DIR}/salmon_index.stdout.log"

    echo "[OK] Salmon index built: ${SALMON_IDX}"
else
    echo "[SKIP] Salmon index already exists."
fi

# ==============================================================================
# 2.  KALLISTO INDEX
# ==============================================================================
if [ ! -f "${KALLISTO_IDX}" ]; then
    echo "=== Building kallisto index ==="
    /usr/bin/time -v \
        -o "${LOG_DIR}/kallisto_index.log" \
        kallisto index \
            -i "${KALLISTO_IDX}" \
            "${TRANSCRIPTOME}" \
        2>&1 | tee "${LOG_DIR}/kallisto_index.stdout.log"
    echo "[OK] kallisto index built: ${KALLISTO_IDX}"
else
    echo "[SKIP] kallisto index already exists."
fi

# ==============================================================================
# 3.  QUANTIFICATION FUNCTION
# ==============================================================================
run_pseudoalignment() {
    local SAMPLE="$1"
    local R1="${RAW_DIR}/${SAMPLE}_1.fastq"
    local R2="${RAW_DIR}/${SAMPLE}_2.fastq"

    if [ ! -f "${R1}" ] || [ ! -f "${R2}" ]; then
        echo "[WARN] FASTQs not found for ${SAMPLE} — skipping."
        return 1
    fi

    # ---- Salmon ----
    local SALMON_OUT="${OUT_DIR}/salmon/${SAMPLE}"
    if [ ! -f "${SALMON_OUT}/quant.sf" ]; then
        echo "[RUN] Salmon — ${SAMPLE}"
        /usr/bin/time -v \
            -o "${LOG_DIR}/salmon_quant_${SAMPLE}.log" \
            salmon quant \
                -i "${SALMON_IDX}" \
                -l A \
                -1 "${R1}" \
                -2 "${R2}" \
                -p "${THREADS}" \
                --gcBias \
                --seqBias \
                --validateMappings \
                -o "${SALMON_OUT}" \
            2>&1 | tee "${LOG_DIR}/salmon_quant_${SAMPLE}.stdout.log"
        echo "[OK] Salmon — ${SAMPLE}"
    else
        echo "[SKIP] Salmon already done for ${SAMPLE}."
    fi

    # ---- kallisto ----
    local KALLISTO_OUT="${OUT_DIR}/kallisto/${SAMPLE}"
    if [ ! -f "${KALLISTO_OUT}/abundance.tsv" ]; then
        echo "[RUN] kallisto — ${SAMPLE}"
        /usr/bin/time -v \
            -o "${LOG_DIR}/kallisto_quant_${SAMPLE}.log" \
            kallisto quant \
                -i "${KALLISTO_IDX}" \
                -o "${KALLISTO_OUT}" \
                -t "${THREADS}" \
                -b "${BOOTSTRAP}" \
                "${R1}" "${R2}" \
            2>&1 | tee "${LOG_DIR}/kallisto_quant_${SAMPLE}.stdout.log"
        echo "[OK] kallisto — ${SAMPLE}"
    else
        echo "[SKIP] kallisto already done for ${SAMPLE}."
    fi
}

# ==============================================================================
# 4.  RUN QUANTIFICATION ON AIM 1 + AIM 2 SAMPLES
# ==============================================================================
AIM1_SAMPLES=( SRR896663 SRR896679 SRR896743 SRR896759 )
AIM2A_SAMPLES=( SRR1643426 SRR1643427 SRR1643434 )
AIM2B_SAMPLES=( )

ALL_QUANT_SAMPLES=( "${AIM1_SAMPLES[@]}" "${AIM2A_SAMPLES[@]}" "${AIM2B_SAMPLES[@]}" )

echo "=== Running Salmon and kallisto (Aim 1 & 2) ==="
for SAMPLE in "${ALL_QUANT_SAMPLES[@]}"; do
    echo "--- ${SAMPLE} ---"
    run_pseudoalignment "${SAMPLE}"
done

echo ""
echo "[DONE] Pseudoalignment complete."
echo "  Salmon results  : ${OUT_DIR}/salmon"
echo "  kallisto results: ${OUT_DIR}/kallisto"
echo "  Logs            : ${LOG_DIR}"
