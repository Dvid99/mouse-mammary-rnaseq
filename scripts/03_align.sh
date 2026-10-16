#!/usr/bin/env bash
# Step 3 — Alignment to GRCm38/mm10 with HISAT2, streamed straight into a
# coordinate-sorted, indexed BAM (avoids writing an ~8 GB intermediate SAM).
source "$(dirname "$0")/utils/common.sh"
require hisat2 samtools

trim_dir="${OUT_DIR}/02_trimmed"
out="${OUT_DIR}/03_alignments"
mkdir -p "$out"

for run in $(run_accessions); do
    log "HISAT2 (-k ${HISAT2_K}): $run"
    hisat2 -p "$THREADS" -k "$HISAT2_K" -x "$HISAT2_INDEX" \
        -U "${trim_dir}/${run}_trimmed.fq.gz" \
        2> "${out}/${run}.hisat2.log" |
        samtools sort -@ "$THREADS" -o "${out}/${run}.sorted.bam" -

    samtools index "${out}/${run}.sorted.bam"
    cat "${out}/${run}.hisat2.log" >&2
done
