#!/usr/bin/env bash
# Step 2 — Adapter and quality trimming with Trim Galore (default settings:
# Phred Q20, auto-detected adapter, min length 20 bp), FastQC on trimmed reads
# and a MultiQC report comparing raw vs trimmed.
source "$(dirname "$0")/utils/common.sh"
require trim_galore fastqc multiqc

trim_dir="${OUT_DIR}/02_trimmed"
qc_dir="${OUT_DIR}/02_qc_trimmed"
mkdir -p "$trim_dir" "$qc_dir"

for run in $(run_accessions); do
    log "Trim Galore: $run"
    trim_galore --cores "$THREADS" --output_dir "$trim_dir" "${RAW_DIR}/${run}.fastq.gz"

    log "FastQC (trimmed): $run"
    fastqc --threads "$THREADS" --outdir "$qc_dir" "${trim_dir}/${run}_trimmed.fq.gz"
done

log "MultiQC report"
multiqc --force --outdir "${OUT_DIR}/multiqc" "${OUT_DIR}/01_qc_raw" "$trim_dir" "$qc_dir"
