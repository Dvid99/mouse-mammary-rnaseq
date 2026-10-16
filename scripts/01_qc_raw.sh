#!/usr/bin/env bash
# Step 1 — Quality control of raw reads with FastQC.
source "$(dirname "$0")/utils/common.sh"
require fastqc

out="${OUT_DIR}/01_qc_raw"
mkdir -p "$out"

for run in $(run_accessions); do
    log "FastQC (raw): $run"
    fastqc --threads "$THREADS" --outdir "$out" "${RAW_DIR}/${run}.fastq.gz"
done
