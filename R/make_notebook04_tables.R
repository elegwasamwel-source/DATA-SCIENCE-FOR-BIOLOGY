# R/make_notebook04_tables.R ---------------------------------------------------
#
# Runs the two notebook 04 steps that save files for notebook 05:
#
#   Exercise 2, step 3 -> data/species_by_outcome.csv  (species ~ outcome,
#                         poscounts size factors)
#   Exercise 2, step 5 -> data/genes_by_outcome.csv    (genes ~ outcome)
#
# Both use the week 18 samples only. Run with:
#   Rscript R/make_notebook04_tables.R

suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(here)
})

source(here("R", "count_tables.R"))
d <- read_counts("instructional")

sample_groups <- read_csv(here("data", "my_sample_groups.csv"),
                          show_col_types = FALSE)
metadata <- d$metadata |> left_join(sample_groups, by = "sample_id")

# Week 18 samples only -- one visit per participant
meta18 <- metadata |> filter(time_point == "week_18")
rownames(meta18) <- meta18$sample_id   # DESeq2 wants the row names to match

# DESeq2 1.50 no longer takes a plain matrix + colData. Build a
# SummarizedExperiment: the count matrix plus its sample table in one object.
make_dds <- function(counts, meta, design_formula) {
  se <- SummarizedExperiment::SummarizedExperiment(
    assays  = S4Vectors::SimpleList(counts = counts),
    colData = S4Vectors::DataFrame(meta)
  )
  DESeqDataSet(se, design = design_formula)
}

# --- Exercise 2, step 3: species ~ outcome, poscounts size factors -----------
# Species present in at least 10% of the week 18 samples
keep_species <- rowSums(d$species[, meta18$sample_id] > 0) >= 0.10 * nrow(meta18)
dds_s <- make_dds(d$species[keep_species, meta18$sample_id], meta18, ~ outcome)

# Spelled out (not DESeq()) so the poscounts size factors survive
dds_s <- estimateSizeFactors(dds_s, type = "poscounts")
dds_s <- estimateDispersions(dds_s)
dds_s <- nbinomWaldTest(dds_s)

# positive log2FoldChange = higher in preterm
res_s <- results(dds_s, contrast = c("outcome", "preterm", "term"))

species_out <- as.data.frame(res_s) |>
  rownames_to_column("feature") |>
  left_join(d$taxonomy, by = c("feature" = "taxon")) |>
  select(feature, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj,
         Genus, Family)

write_csv(species_out, here("data", "species_by_outcome.csv"))

# --- Exercise 2, step 5: genes ~ outcome -------------------------------------
# Genes with at least 10 counts in at least a quarter of the week 18 samples
keep_genes <- rowSums(d$genes[, meta18$sample_id] >= 10) >= nrow(meta18) / 4
dds_g <- make_dds(d$genes[keep_genes, meta18$sample_id], meta18, ~ outcome)
dds_g <- DESeq(dds_g)

# positive log2FoldChange = higher in preterm
res_g <- results(dds_g, contrast = c("outcome", "preterm", "term"))

genes_out <- as.data.frame(res_g) |>
  rownames_to_column("feature") |>
  left_join(d$gene_info, by = c("feature" = "ensembl_gene_id")) |>
  select(feature, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj,
         symbol, gene_name)

write_csv(genes_out, here("data", "genes_by_outcome.csv"))

# --- Numbers to check the notebook's questions against -----------------------
cat("\n===== species_by_outcome.csv =====\n")
cat("Species tested:          ", nrow(species_out), "\n")
cat("Significant species:     ",
    sum(species_out$padj < 0.05, na.rm = TRUE), "\n")

cat("\n===== genes_by_outcome.csv =====\n")
cat("Genes tested:            ", nrow(genes_out), "\n")
cat("Significant genes:       ",
    sum(genes_out$padj < 0.05, na.rm = TRUE), "\n")
cat("  higher in preterm:     ",
    sum(genes_out$padj < 0.05 & genes_out$log2FoldChange > 0, na.rm = TRUE), "\n")
cat("  higher in term:        ",
    sum(genes_out$padj < 0.05 & genes_out$log2FoldChange < 0, na.rm = TRUE), "\n")
cat("More than two-fold:      ",
    sum(genes_out$padj < 0.05 & abs(genes_out$log2FoldChange) > 1, na.rm = TRUE), "\n")
