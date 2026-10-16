#!/usr/bin/env bash
# Step 4 — Alignment statistics with SAMtools: flagstat plus the distribution
# of SAM FLAG and MAPQ values (useful to see unique vs multi-mapped reads).
source "$(dirname "$0")/utils/common.sh"
require samtools

aln="${OUT_DIR}/03_alignments"
out="${OUT_DIR}/04_alignment_stats"
mkdir -p "$out"

for run in $(run_accessions); do
    bam="${aln}/${run}.sorted.bam"
    log "SAMtools stats: $run"

    samtools flagstat -@ "$THREADS" "$bam" > "${out}/${run}.flagstat.txt"

    # count  FLAG
    samtools view "$bam" | cut -f2 | sort | uniq -c | sort -rn > "${out}/${run}.flag_counts.txt"
    # count  MAPQ  (HISAT2: 60 = unique, 1 = multi-mapped, 0 = unmapped/ambiguous)
    samtools view "$bam" | cut -f5 | sort -n | uniq -c > "${out}/${run}.mapq_counts.txt"
done
