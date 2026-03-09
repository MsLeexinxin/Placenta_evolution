#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(dplyr)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Error: At least 3 arguments are required: <input.rds> <output_dir> <group_map.csv> <annotation_column>")
}

rds_file   <- args[1]
output_dir <- args[2]
map_file   <- args[3]
target_col <- args[4]

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

seurat_obj <- readRDS(rds_file)

# Extract normalized expression matrix 
counts <- LayerData(seurat_obj, assay = "RNA", layer = "data") %>%
  as.data.frame() %>%
  rownames_to_column(var = "Gene")

# Extract basic metadata
if (!target_col %in% colnames(seurat_obj@meta.data)) {
  stop(sprintf("Error: Column '%s' not found in Seurat object!", target_col))
}

metadata <- data.frame(
  Cell = colnames(seurat_obj),
  cell_type = as.character(seurat_obj@meta.data[[target_col]]),
  stringsAsFactors = FALSE
)

# Add Fetal/Maternal group info (if mapping file is provided)
if (!is.na(map_file) && file.exists(map_file)) {
  cat(sprintf("[3/3] Applying Fetal/Maternal group mapping (%s)...\n", map_file))
  group_map <- read.csv(map_file, stringsAsFactors = FALSE, strip.white = TRUE)
  
  if (!all(c("cell_type", "group") %in% colnames(group_map))) {
    stop("Error: Mapping file must contain 'cell_type' and 'group' columns!")
  }
  
  # Merge group information
  metadata <- metadata %>%
    left_join(group_map, by = "cell_type")
  
  # Check for unmapped cell types
  unmapped <- unique(metadata$cell_type[is.na(metadata$group)])
  if (length(unmapped) > 0) {
    warning(sprintf("\nWarning: The following cell types were not found in the mapping file and are marked as NA:\n  %s", paste(unmapped, collapse=", ")))
  }
} else {
  cat("mapping.csv not provided or not found. Only Cell and cell_type columns will be output.\n")
}

# Save files
output_prefix <- tools::file_path_sans_ext(basename(rds_file))
output_counts_file <- file.path(output_dir, paste0(output_prefix, "_counts.txt"))
output_meta_file   <- file.path(output_dir, paste0(output_prefix, "_meta.txt"))
fwrite(counts, file = output_counts_file, row.names = FALSE, sep = "\t")
fwrite(metadata, file = output_meta_file, row.names = FALSE, sep = "\t")
