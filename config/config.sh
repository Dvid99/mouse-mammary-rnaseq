# Pipeline configuration — edit these paths for your system.
# All paths are relative to the repository root (scripts cd there automatically).

SAMPLES="config/samples.tsv"        # sample sheet (run_accession column = FASTQ name)
RAW_DIR="data/raw"                  # {run_accession}.fastq.gz
REF_DIR="data/reference"

HISAT2_INDEX="${REF_DIR}/mm10/genome"                 # HISAT2 index prefix (GRCm38/mm10)
GTF="${REF_DIR}/gencode.vM10.annotation.gtf"          # GENCODE M10 annotation

OUT_DIR="results"
THREADS=4

HISAT2_K=1             # max alignments per read (-k); 1 is appropriate for HTSeq-count
HTSEQ_STRANDED="no"    # TruSeq RNA v2 (non-stranded) -> "no"
QUALIMAP_MEM="4000M"

# Optional: run only a subset of samples, e.g. SAMPLE_FILTER="SRR1552444"
SAMPLE_FILTER="${SAMPLE_FILTER:-}"
