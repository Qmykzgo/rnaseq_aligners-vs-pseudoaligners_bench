#!/bin/bash
# ==============================================================================
# 00_setup_conda_envs.sh
# Purpose : Create reproducible Conda environments for the RNA-seq benchmark.
# Environments:
#   rnaseq_bench_env   — aligners, quantifiers, sra-tools, samtools, gffread
#   rnaseq_metrics_env — Python data-science stack for metric aggregation
#
# Usage: bash 00_setup_conda_envs.sh
# ==============================================================================
set -euo pipefail

# ---------- Locate and source Conda ----------
CONDA_BASE="$(conda info --base 2>/dev/null || echo "${HOME}/miniconda3")"
source "${CONDA_BASE}/etc/profile.d/conda.sh"

echo "============================================================"
echo "  RNA-seq Benchmark — Conda Environment Setup"
echo "  Conda base : ${CONDA_BASE}"
echo "============================================================"

# ---------- Helper: skip if env already exists ----------
env_exists() { conda env list | awk '{print $1}' | grep -qx "$1"; }

# ==============================================================================
# Environment 1: Aligners, quantifiers, data-download utilities
# ==============================================================================
ENV1="rnaseq_bench_env"
if env_exists "${ENV1}"; then
    echo "[SKIP] Environment '${ENV1}' already exists."
else
    echo "[INFO] Creating '${ENV1}' …"
    conda create -n "${ENV1}" -y \
        -c conda-forge \
        -c bioconda \
        python=3.10 \
        star=2.7.10b \
        hisat2=2.2.1 \
        subread=2.0.6 \
        salmon=1.10.3 \
        kallisto=0.50.1 \
        sra-tools=3.1.1 \
        samtools=1.20 \
        htslib=1.20 \
        gffread=0.12.7 \
        wget \
        unzip
    echo "[OK] '${ENV1}' created."
fi

# ==============================================================================
# Environment 2: Python data-science for downstream metric parsing
# ==============================================================================
ENV2="rnaseq_metrics_env"
if env_exists "${ENV2}"; then
    echo "[SKIP] Environment '${ENV2}' already exists."
else
    echo "[INFO] Creating '${ENV2}' …"
    conda create -n "${ENV2}" -y \
        -c conda-forge \
        python=3.10 \
        pandas=2.1.4 \
        numpy=1.26.3 \
        scipy=1.12.0 \
        matplotlib=3.8.2 \
        seaborn=0.13.2 \
        statsmodels=0.14.1
    echo "[OK] '${ENV2}' created."
fi

echo ""
echo "[DONE] All environments ready."
printf "  %-30s %s\n" "Activate aligners / tools :" "conda activate ${ENV1}"
printf "  %-30s %s\n" "Activate metrics / Python  :" "conda activate ${ENV2}"
