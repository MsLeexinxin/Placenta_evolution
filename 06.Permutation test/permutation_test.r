#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(dplyr)
  library(tools)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Error: At least 5 arguments are required.
Usage: Rscript 08_permutation_test.R <hemo_genes.txt> <bg_genes.txt> <cancer_dir> <out_summary.csv> <out_intersections.csv> [iterations] [seed]")
}

HEMO_FILE     <- args[1]
BG_FILE       <- args[2]
CANCER_DIR    <- args[3]
OUT_SUMMARY   <- args[4]
OUT_INTERSECT <- args[5]
ITERATIONS <- ifelse(length(args) >= 6, as.numeric(args[6]), 10000)
SET_SEED   <- ifelse(length(args) >= 7, as.numeric(args[7]), 42)

set.seed(SET_SEED)

# Loading background genes
if (!file.exists(HEMO_FILE) || !file.exists(BG_FILE)) {
  stop("Error: Hemochorial or background genes file not found!")
}

hemo_up_genes <- unique(as.character(read.table(HEMO_FILE, header = FALSE)$V1))
all_background_genes <- unique(as.character(read.table(BG_FILE, header = FALSE)$V1))

# Loadign cancer up-regulated genes
cancer_files <- list.files(CANCER_DIR, pattern = "\\.txt$", full.names = TRUE)
if (length(cancer_files) == 0) {
  stop(sprintf("Error: No .txt files found in directory '%s'", CANCER_DIR))
}

cancer_names <- file_path_sans_ext(basename(cancer_files))

summary_results <- list()
intersection_genes_list <- list()

# Porcessing each cancer file
for (i in seq_along(cancer_files)) {
  current_cancer <- cancer_names[i]  
  disease_genes <- unique(as.character(read.table(cancer_files[i], header = FALSE)$V1))
  overlap_genes <- intersect(disease_genes, hemo_up_genes)
  observed_count <- length(overlap_genes)
  is_in_disease <- all_background_genes %in% disease_genes
  
  extreme_count <- 0
  for (j in 1:ITERATIONS) {
    random_indices <- sample(seq_along(all_background_genes), size = length(hemo_up_genes), replace = FALSE)
    random_overlap_count <- sum(is_in_disease[random_indices])
    
    if (random_overlap_count >= observed_count) {
      extreme_count <- extreme_count + 1
    }
  }

  p_val <- (extreme_count + 1) / (ITERATIONS + 1)
  
  summary_results[[current_cancer]] <- data.frame(
    `Cancer type` = current_cancer,
    `Number of up-regulated genes in cancer cells` = length(disease_genes),
    `Number of up-regulated genes in hemochorial cells` = length(hemo_up_genes),
    `Number of overlapped genes` = observed_count,
    `Number of more than overlapping sweep` = extreme_count,
    `P value` = p_val,
    check.names = FALSE
  ) 
  intersection_genes_list[[current_cancer]] <- overlap_genes
}

# Save results
# Calculate padj (BH method)
final_summary <- do.call(rbind, summary_results)
final_summary$padj <- p.adjust(final_summary$`P value`, method = "BH")
write.csv(final_summary, OUT_SUMMARY, row.names = FALSE)

# Write the list of intersecting genes
max_len <- max(sapply(intersection_genes_list, length))
padded_list <- lapply(intersection_genes_list, function(x) {
  if (length(x) == 0) return(rep("", max_len))
  c(x, rep("", max_len - length(x)))
})
gene_table <- as.data.frame(do.call(cbind, padded_list))
write.csv(gene_table, OUT_INTERSECT, row.names = FALSE)
cat(sprintf("Intersection genes list saved to: %s\n", OUT_INTERSECT))
