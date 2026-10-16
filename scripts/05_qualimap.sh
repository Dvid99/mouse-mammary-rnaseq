#!/usr/bin/env bash
# Step 5 — RNA-seq specific alignment QC with QualiMap (genomic origin of
# reads, 5'/3' coverage bias, junctions).
source "$(dirname "$0")/utils/common.sh"
require qualimap

aln="${OUT_DIR}/03_alignments"
out="${OUT_DIR}/05_qualimap"

for run in $(run_accessions); do
    log "QualiMap rnaseq: $run"
    qualimap rnaseq \
        -bam "${aln}/${run}.sorted.bam" \
        -gtf "$GTF" \
        --java-mem-size="$QUALIMAP_MEM" \
        -outdir "${out}/${run}"
done
