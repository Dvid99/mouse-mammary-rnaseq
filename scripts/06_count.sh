#!/usr/bin/env bash
# Step 6 — Gene-level read counting with HTSeq-count (union mode, unstranded),
# then a summary plot of assigned vs. unassigned reads.
source "$(dirname "$0")/utils/common.sh"
require htseq-count python

aln="${OUT_DIR}/03_alignments"
out="${OUT_DIR}/06_counts"
mkdir -p "$out"

for run in $(run_accessions); do
    log "HTSeq-count: $run"
    htseq-count \
        --format bam \
        --order pos \
        --stranded "$HTSEQ_STRANDED" \
        --idattr gene_id \
        "${aln}/${run}.sorted.bam" "$GTF" \
        > "${out}/${run}.counts.tsv"
done

log "Plotting HTSeq-count assignment summary"
python scripts/07_plot_count_summary.py "${out}"/*.counts.tsv \
    --out "${out}/htseq_assignment_summary.png" \
    --table "${out}/htseq_assignment_summary.tsv"
