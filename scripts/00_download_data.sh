#!/usr/bin/env bash
# Step 0 — Download public input data (nothing is stored in the repository).
#   bash scripts/00_download_data.sh            # all 12 FASTQs (~1.5 GB each) + reference + counts
#   SAMPLE_FILTER=SRR1552444 bash scripts/00_download_data.sh   # one sample only
#   bash scripts/00_download_data.sh --test     # 1,000-read test FASTQs instead of full reads
#                                               # (Galaxy Training, Zenodo) + reference + counts
source "$(dirname "$0")/utils/common.sh"
require wget

mkdir -p "$RAW_DIR" "$REF_DIR" data/test

if [[ "${1:-}" == "--test" ]]; then
    for run in $(run_accessions); do
        log "Test FASTQ: $run"
        wget -nc -q -O "data/test/${run}.fastq.gz" \
            "https://zenodo.org/records/4249516/files/${run}.fastq.gz"
    done
else
    # Raw reads from ENA (Illumina HiSeq 2000, 100 bp single-end)
    for run in $(run_accessions); do
        log "FASTQ: $run"
        wget -nc -P "$RAW_DIR" \
            "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${run:0:6}/00${run: -1}/${run}/${run}.fastq.gz"
    done
fi

# GENCODE M10 annotation (GRCm38)
if [[ ! -s "$GTF" ]]; then
    log "GTF: GENCODE vM10"
    wget -O - "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M10/gencode.vM10.annotation.gtf.gz" |
        gunzip > "$GTF"
fi

# Pre-built HISAT2 index for GRCm38/mm10 (~3.5 GB)
if [[ ! -s "${HISAT2_INDEX}.1.ht2" ]]; then
    log "HISAT2 index: mm10"
    wget -O - "https://genome-idx.s3.amazonaws.com/hisat/mm10_genome.tar.gz" | tar -xz -C "$REF_DIR"
fi

# Gene-wise count matrix used by the downstream (edgeR) analysis
if [[ ! -s data/GSE60450_Lactation-GenewiseCounts.txt ]]; then
    log "GEO count matrix: GSE60450"
    wget -O - "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE60450&format=file&file=GSE60450_Lactation-GenewiseCounts.txt.gz" |
        gunzip > data/GSE60450_Lactation-GenewiseCounts.txt
fi
