#!/usr/bin/env Rscript
# Step 8 — Differential expression with edgeR (quasi-likelihood GLM).
#
# Input : GSE60450 gene-wise count matrix (12 samples, Entrez gene IDs)
# Output: QC plots, MDS, BCV, DE tables, volcano plot, Venn diagram and a
#         heatmap of genes up-regulated in lactation in BOTH cell types.
#
# Usage : Rscript scripts/08_differential_expression.R

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(VennDiagram)
})

# ---- Paths & parameters ------------------------------------------------------
counts_file <- "data/GSE60450_Lactation-GenewiseCounts.txt"
samples_file <- "config/samples.tsv"
out_dir <- "results/08_differential_expression"
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

FDR_CUTOFF <- 0.05
LOGFC_CUTOFF <- 1

save_png <- function(name, expr, width = 1800, height = 1200, res = 200) {
  png(file.path(fig_dir, name), width = width, height = height, res = res)
  on.exit(dev.off())
  force(expr)
}

# ---- 1. Read counts and sample sheet -----------------------------------------
raw <- read.delim(counts_file, row.names = 1, check.names = FALSE)
counts <- as.matrix(raw[, -1])                  # drop the 'Length' column

samples <- read.delim(samples_file)
# GEO column names look like "MCL1-DG_BC2CTUACXX_ACTTGA_L002_R1" -> "MCL1.DG"
geo_names <- gsub("-", ".", sub("_.*$", "", colnames(counts)), fixed = TRUE)
samples <- samples[match(geo_names, samples$geo_sample_name), ]
if (anyNA(samples$sample_id)) {
  stop("Count-matrix columns not found in ", samples_file, ": ",
       paste(colnames(counts)[is.na(samples$sample_id)], collapse = ", "))
}
colnames(counts) <- samples$sample_id
group <- factor(paste(samples$cell_type, samples$status, sep = "."))

message("Genes in raw matrix: ", nrow(counts))

# Gene with the highest count in each sample
max_counts <- apply(counts, 2, max)
symbols <- mapIds(org.Mm.eg.db, keys = rownames(counts), column = "SYMBOL",
                  keytype = "ENTREZID", multiVals = "first")
message("Entrez IDs without SYMBOL: ", sum(is.na(symbols)))
max_gene <- symbols[apply(counts, 2, which.max)]

df_max <- data.frame(sample = factor(names(max_counts), levels = names(max_counts)),
                     max = max_counts, gene = max_gene)
p <- ggplot(df_max, aes(sample, max, fill = gene)) +
  geom_col() +
  scale_y_continuous(labels = scales::label_comma()) +
  labs(title = "Highest gene count per sample", x = NULL, y = "Counts", fill = "Gene") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(fig_dir, "max_counts_per_sample.png"), p, width = 8, height = 5, dpi = 200)

pct_detected <- round(colMeans(counts > 0) * 100, 2)
write.table(data.frame(sample = names(pct_detected), pct_genes_detected = pct_detected),
            file.path(out_dir, "pct_genes_detected.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE)

# ---- 2. DGEList and low-expression filter ------------------------------------
genes <- data.frame(ENTREZID = rownames(counts), SYMBOL = symbols,
                    Length = raw$Length, row.names = rownames(counts))
dge <- DGEList(counts = counts, group = group, genes = genes)

keep <- filterByExpr(dge)
message(sprintf("Genes kept by filterByExpr: %d / %d (%.2f%%)",
                sum(keep), length(keep), mean(keep) * 100))
dge <- dge[keep, , keep.lib.sizes = FALSE]

save_png("logcpm_before_tmm.png",
         boxplot(cpm(dge, log = TRUE), las = 2, col = "steelblue",
                 main = "logCPM before TMM normalisation", ylab = "log2 CPM",
                 cex.axis = 0.7))

# ---- 3. TMM normalisation and MDS --------------------------------------------
dge <- calcNormFactors(dge, method = "TMM")
write.table(dge$samples, file.path(out_dir, "tmm_norm_factors.tsv"), sep = "\t",
            quote = FALSE, col.names = NA)

logcpm <- cpm(dge, log = TRUE)
save_png("logcpm_after_tmm.png",
         boxplot(logcpm, las = 2, col = "steelblue",
                 main = "logCPM after TMM normalisation", ylab = "log2 CPM",
                 cex.axis = 0.7))

status_col <- c(virgin = "blue", pregnant = "red", lactate = "green3")
celltype_pch <- c(basal = 1, luminal = 17)
save_png("mds_plot.png", width = 2000, expr = {
  par(mar = c(5, 5, 4, 9), xpd = TRUE)
  plotMDS(dge, col = status_col[samples$status], pch = celltype_pch[samples$cell_type],
          cex = 1.4, main = "MDS plot")
  legend("topright", inset = c(-0.22, 0), bty = "n", cex = 0.9,
         legend = c(names(status_col), names(celltype_pch)),
         col = c(status_col, "black", "black"),
         pch = c(16, 16, 16, celltype_pch))
})

# ---- 4. Dispersion estimation and QL GLM fit ---------------------------------
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

dge <- estimateDisp(dge, design)
message("Common dispersion: ", signif(dge$common.dispersion, 4),
        " (BCV = ", signif(sqrt(dge$common.dispersion), 3), ")")
save_png("bcv_plot.png", plotBCV(dge))

# edgeR >= 4.0 changed the default quasi-likelihood method; legacy = TRUE keeps
# the classic method (edgeR 3.x) used for the results reported in the README.
fit <- if (packageVersion("edgeR") >= "4.0.0") {
  glmQLFit(dge, design, legacy = TRUE)
} else {
  glmQLFit(dge, design)
}

# ---- 5. Contrasts: lactation vs pregnancy in each cell type ------------------
run_contrast <- function(expr, label) {
  con <- makeContrasts(contrasts = expr, levels = design)
  res <- topTags(glmQLFTest(fit, contrast = con), n = Inf)$table
  write.table(res, file.path(out_dir, paste0("DE_", label, ".tsv")), sep = "\t",
              quote = FALSE, row.names = FALSE)
  degs <- res[res$FDR <= FDR_CUTOFF & abs(res$logFC) >= LOGFC_CUTOFF, ]
  message(sprintf("%s: %d DEGs (%d up, %d down)", label, nrow(degs),
                  sum(degs$logFC > 0), sum(degs$logFC < 0)))
  list(table = res,
       up = degs$ENTREZID[degs$logFC > 0],
       down = degs$ENTREZID[degs$logFC < 0])
}

luminal <- run_contrast("luminal.lactate - luminal.pregnant", "luminal_lactate_vs_pregnant")
basal <- run_contrast("basal.lactate - basal.pregnant", "basal_lactate_vs_pregnant")

# Volcano plot (luminal), labelling the 3 genes with the largest |logFC|
res <- luminal$table
res$status <- "Not significant"
res$status[res$FDR <= FDR_CUTOFF & res$logFC >= LOGFC_CUTOFF] <- "Up"
res$status[res$FDR <= FDR_CUTOFF & res$logFC <= -LOGFC_CUTOFF] <- "Down"
top3 <- head(res[order(abs(res$logFC), decreasing = TRUE), ], 3)

p <- ggplot(res, aes(logFC, -log10(FDR), colour = status)) +
  geom_point(alpha = 0.5, size = 1) +
  scale_colour_manual(values = c(Up = "red", Down = "blue", `Not significant` = "grey70")) +
  geom_vline(xintercept = c(-LOGFC_CUTOFF, LOGFC_CUTOFF), linetype = "dashed") +
  geom_hline(yintercept = -log10(FDR_CUTOFF), linetype = "dashed") +
  geom_label_repel(data = top3, aes(label = SYMBOL), size = 3, show.legend = FALSE) +
  labs(title = "Luminal cells: lactation vs pregnancy",
       x = "log2 fold change", y = "-log10(FDR)", colour = NULL) +
  theme_classic()
ggsave(file.path(fig_dir, "volcano_luminal_lactate_vs_pregnant.png"), p,
       width = 7, height = 5.5, dpi = 200)

# Venn diagram: overlap of up/down DEGs between cell types
invisible(futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger"))
invisible(venn.diagram(
  x = list(`Luminal UP` = luminal$up, `Luminal DOWN` = luminal$down,
           `Basal UP` = basal$up, `Basal DOWN` = basal$down),
  filename = file.path(fig_dir, "venn_degs_luminal_vs_basal.png"),
  imagetype = "png", col = c("red", "blue", "orange", "green"),
  fill = c("red", "blue", "orange", "green"), alpha = 0.3, cex = 1.2,
  fontface = "bold", main = "Lactation vs pregnancy DEGs"
))

# ---- 6. Heatmap of genes up-regulated in BOTH cell types ---------------------
shared_up <- intersect(luminal$up, basal$up)
message("Genes up-regulated in lactation in both cell types: ", length(shared_up))
writeLines(shared_up, file.path(out_dir, "shared_up_entrez.txt"))

if (length(shared_up) >= 2) {
  mat <- logcpm[shared_up, , drop = FALSE]
  rownames(mat) <- ifelse(is.na(dge$genes[shared_up, "SYMBOL"]), shared_up,
                          dge$genes[shared_up, "SYMBOL"])
  annotation_col <- data.frame(CellType = samples$cell_type, Status = samples$status,
                               row.names = samples$sample_id)
  pheatmap(mat, scale = "row", show_rownames = FALSE,
           clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",
           clustering_method = "complete", annotation_col = annotation_col,
           main = "Genes up-regulated in lactation (basal & luminal)", fontsize_col = 8,
           filename = file.path(fig_dir, "heatmap_shared_upregulated.png"),
           width = 8, height = 7)
} else {
  message("Fewer than 2 shared up-regulated genes: heatmap skipped")
}

# Objects needed by the enrichment step
saveRDS(list(universe = rownames(dge), shared_up = shared_up),
        file.path(out_dir, "enrichment_input.rds"))

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
message("Done. Results in ", out_dir)
