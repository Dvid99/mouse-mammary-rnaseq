#!/usr/bin/env bash
# Exploratory summary of a GTF annotation with standard Unix tools:
# annotation sources, feature types, genes per chromosome and longest genes.
# Usage: bash scripts/utils/gtf_summary.sh [annotation.gtf]
source "$(dirname "$0")/common.sh"

gtf="${1:-$GTF}"
out="${OUT_DIR}/gtf_summary"
mkdir -p "$out"

body() { grep -v '^#' "$gtf"; }

log "Annotation sources (column 2)"
body | cut -f2 | sort | uniq -c | sort -rn | tee "${out}/sources.txt"

log "Feature types (column 3)"
body | cut -f3 | sort | uniq -c | sort -rn | tee "${out}/features.txt"

log "Genes per chromosome"
body | awk -F'\t' '$3 == "gene" { print $1 }' | sort | uniq -c | sort -rn \
    | tee "${out}/genes_per_chromosome.txt"

log "Top 5 longest genes (length, chr, start, end, gene_name, gene_id)"
body | awk -F'\t' '$3 == "gene" {
        match($9, /gene_name "[^"]+"/); name = substr($9, RSTART + 11, RLENGTH - 12)
        match($9, /gene_id "[^"]+"/);   id   = substr($9, RSTART + 9,  RLENGTH - 10)
        print $5 - $4 + 1, $1, $4, $5, name, id
    }' OFS='\t' | sort -k1,1nr | head -5 | tee "${out}/longest_genes.txt"
