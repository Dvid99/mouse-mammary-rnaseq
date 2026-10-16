#!/usr/bin/env Rscript
# Step 9 — GO enrichment (Biological Process and Molecular Function) with topGO
# for the genes up-regulated in lactation in both basal and luminal cells.
# Universe = all genes that passed filterByExpr in step 8.
#
# Usage : Rscript scripts/09_go_enrichment.R

suppressPackageStartupMessages({
  library(topGO)
  library(org.Mm.eg.db)
})

in_dir <- "results/08_differential_expression"
out_dir <- "results/09_go_enrichment"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

input <- readRDS(file.path(in_dir, "enrichment_input.rds"))
gene_list <- factor(as.integer(input$universe %in% input$shared_up), levels = c(0, 1))
names(gene_list) <- input$universe
message("Genes of interest: ", sum(gene_list == 1), " / universe: ", length(gene_list))

run_go <- function(ontology, top_n = 10) {
  go_data <- new("topGOdata", ontology = ontology, allGenes = gene_list,
                 annot = annFUN.org, mapping = "org.Mm.eg.db", ID = "entrez")
  result <- runTest(go_data, algorithm = "classic", statistic = "fisher")
  tab <- GenTable(go_data, classicFisher = result, topNodes = top_n, numChar = 80)
  write.table(tab, file.path(out_dir, paste0("topGO_", ontology, "_top", top_n, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  print(tab)
  invisible(tab)
}

run_go("BP")
run_go("MF")

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
