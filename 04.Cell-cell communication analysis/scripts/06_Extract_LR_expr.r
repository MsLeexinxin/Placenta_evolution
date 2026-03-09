#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Error: Exactly 4 arguments required: <data_dir> <lr_database> <output_dir> <species_list>")
}

DATA_DIR <- args[1]
LR_FILE  <- args[2]
OUTPUT_DIR <- args[3]
SPECIES_LIST <- strsplit(args[4], ",")[[1]]

if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

lr_pairs <- read.csv(LR_FILE)
ligand_genes <- unique(lr_pairs$Ligand)
receptor_genes <- unique(lr_pairs$Receptor)

for (species in SPECIES_LIST) {
  cat(paste("Processing species:", species, "\n"))
  
  counts_path <- file.path(DATA_DIR, paste0(species, "_counts.rds"))
  metadata_path <- file.path(DATA_DIR, paste0(species, "_meta.txt"))
  
  if (!file.exists(counts_path) || !file.exists(metadata_path)) {
      cat(paste("Warning: Missing data for", species, "Skipping.\n"))
      next
  }
  
  counts_matrix <- readRDS(counts_path)
  metadata <- read.csv(metadata_path, row.names = 1)
  
  if (!all(colnames(counts_matrix) %in% rownames(metadata))) {
    cat(paste("Error: Cell IDs mismatch for", species, "! Skipping.\n"))
    next
  }
  
  # Calculate CPM & zTPM stats
  cpm_matrix <- Matrix::t(Matrix::t(counts_matrix) / Matrix::colSums(counts_matrix)) * 1e6
  log2_cpm_matrix <- log2(cpm_matrix + 1)
  gene_means_log2_cpm <- Matrix::rowMeans(log2_cpm_matrix)
  
  stats_df <- data.frame(
      species = species, 
      mu = mean(gene_means_log2_cpm, na.rm = TRUE), 
      sigma = sd(gene_means_log2_cpm, na.rm = TRUE)
  )
  write.csv(stats_df, file.path(OUTPUT_DIR, paste0(species, "_tpm_stats.csv")), row.names = FALSE)
  
  fetal_cells <- rownames(metadata[metadata$group == "Fetal", ])
  maternal_cells <- rownames(metadata[metadata$group == "Maternal", ])
  
  # Process Ligands (Fetal)
  existing_ligands <- intersect(ligand_genes, rownames(cpm_matrix))
  if (length(fetal_cells) == 0 || length(existing_ligands) == 0) {
    cat("Warning: No Fetal cells or matched ligands.\n")
    ligand_output <- merge(lr_pairs, data.frame(Ligand=character(), TPM=numeric()), by = "Ligand", all.x = TRUE)
  } else {
    fetal_cpm <- cpm_matrix[existing_ligands, fetal_cells, drop = FALSE]
    fetal_metadata <- metadata[fetal_cells, ]
    avg_fetal_expr <- sapply(unique(fetal_metadata$cell_type), function(ct) {
      Matrix::rowMeans(fetal_cpm[, rownames(fetal_metadata[fetal_metadata$cell_type == ct, ]), drop = FALSE])
    })
    ligand_repr_expr <- if (is.null(dim(avg_fetal_expr))) avg_fetal_expr else apply(avg_fetal_expr, 1, max, na.rm = TRUE)
    ligand_output <- merge(lr_pairs, data.frame(Ligand = names(ligand_repr_expr), TPM = ligand_repr_expr), by = "Ligand", all.x = TRUE)
  }
  write.csv(ligand_output, file.path(OUTPUT_DIR, paste0(species, "_ligand_expression.csv")), row.names = FALSE)
  
  # Process Receptors (Maternal)
  existing_receptors <- intersect(receptor_genes, rownames(cpm_matrix))
  if (length(maternal_cells) == 0 || length(existing_receptors) == 0) {
    cat("Warning: No Maternal cells or matched receptors.\n")
    receptor_output <- merge(lr_pairs, data.frame(Receptor=character(), TPM=numeric()), by = "Receptor", all.x = TRUE)
  } else {
    maternal_cpm <- cpm_matrix[existing_receptors, maternal_cells, drop = FALSE]
    maternal_metadata <- metadata[maternal_cells, ]
    avg_maternal_expr <- sapply(unique(maternal_metadata$cell_type), function(ct) {
      Matrix::rowMeans(maternal_cpm[, rownames(maternal_metadata[maternal_metadata$cell_type == ct, ]), drop = FALSE])
    })
    receptor_repr_expr <- if (is.null(dim(avg_maternal_expr))) avg_maternal_expr else apply(avg_maternal_expr, 1, max, na.rm = TRUE)
    receptor_output <- merge(lr_pairs, data.frame(Receptor = names(receptor_repr_expr), TPM = receptor_repr_expr), by = "Receptor", all.x = TRUE)
  }
  write.csv(receptor_output, file.path(OUTPUT_DIR, paste0(species, "_receptor_expression.csv")), row.names = FALSE)
}