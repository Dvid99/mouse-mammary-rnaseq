#!/usr/bin/env bash
# Shared helpers sourced by every upstream script.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"
# shellcheck source=../../config/config.sh
source config/config.sh

log() { printf '[%s] %s\n' "$(date '+%H:%M:%S')" "$*" >&2; }

# Prints the run accessions to process (one per line), honouring SAMPLE_FILTER.
run_accessions() {
    awk -F'\t' 'NR > 1 && $2 != "" { print $2 }' "$SAMPLES" |
        { if [[ -n "$SAMPLE_FILTER" ]]; then grep -Fx -e "$SAMPLE_FILTER"; else cat; fi; }
}

require() {
    for tool in "$@"; do
        command -v "$tool" >/dev/null 2>&1 || { log "ERROR: '$tool' not found. Did you 'conda activate rnaseq-upstream'?"; exit 1; }
    done
}
