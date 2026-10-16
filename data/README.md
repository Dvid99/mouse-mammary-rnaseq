# Data

No data is stored in this repository. Everything is public and can be fetched with
`bash scripts/00_download_data.sh`.

| Input | Source | Used by |
|---|---|---|
| Raw reads, 12 samples (100 bp single-end, Illumina HiSeq 2000) | SRA/ENA `SRR1552444`–`SRR1552455` ([GSE60450](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60450)) | steps 1–7 |
| HISAT2 index, GRCm38/mm10 | [HISAT2 pre-built indexes](https://daehwankimlab.github.io/hisat2/download/) | step 3 |
| Annotation, GENCODE M10 | [GENCODE](https://www.gencodegenes.org/mouse/release_M10.html) | steps 5–6 |
| Gene-wise count matrix (Entrez IDs, featureCounts) | GEO supplementary file `GSE60450_Lactation-GenewiseCounts.txt.gz` | steps 8–9 |
| Test FASTQs (1,000 reads each) | [Galaxy Training Network, Zenodo 4249516](https://zenodo.org/records/4249516) | `-profile test` |

Sample metadata (cell type, physiological status, GEO name) is in
[`config/samples.tsv`](../config/samples.tsv).

Expected layout after downloading:

```
data/
├── raw/SRR15524{44..55}.fastq.gz
├── reference/
│   ├── mm10/genome.*.ht2
│   └── gencode.vM10.annotation.gtf
└── GSE60450_Lactation-GenewiseCounts.txt
```
