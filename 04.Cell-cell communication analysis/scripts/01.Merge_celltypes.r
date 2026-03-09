#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Error: 2 arguments are required: <config.csv> <Output_merged.rds>")
}

config_file <- args[1]
output_path <- args[2]

if (!file.exists(config_file)) {
  stop(paste("Configuration file not found:", config_file))
}

# Read the condig table
config <- read.csv(config_file, stringsAsFactors = FALSE, strip.white = TRUE)
objects_to_merge <- list()

for (i in seq_len(nrow(config))) {
  obj_name <- config$object_name[i]
  file_path <- config$file_path[i]
  meta_col <- config$metadata_col[i]
  keep_clusters <- config$keep_clusters[i]
    
  if (!file.exists(file_path)) {
    warning(paste("File not found, skipping:", file_path))
    next
  }
  
  # loading seurat data
  seurat_data <- readRDS(file_path)
  
  # Check metadata column
  if (!meta_col %in% colnames(seurat_data@meta.data)) {
    warning(sprintf("Metadata column '%s' not found in %s! Skipping.", meta_col, obj_name))
    next
  }  
  # extract data
  if (tolower(keep_clusters) != "all" && keep_clusters != "") {  
    target_clusters <- unlist(strsplit(keep_clusters, "\\|"))
    cat(sprintf("Subsetting for clusters: %s\n", paste(target_clusters, collapse = ", ")))   
    Idents(seurat_data) <- meta_col
    seurat_data <- subset(seurat_data, idents = target_clusters)
  } else {
    cat("Keeping all cells.\n")
  }
  
  # Create columns for CellPhoneDB
  seurat_data$cell_type_cpdb <- seurat_data@meta.data[[meta_col]]
  objects_to_merge[[obj_name]] <- seurat_data
}

# Merging processed objects
if (length(objects_to_merge) == 0) {
  stop("Error: No valid objects were successfully processed. Check your config file and paths.")
} else if (length(objects_to_merge) == 1) {
  merged_seurat <- objects_to_merge[[1]]
  cat("Only one object processed. No merge needed.\n")
} else {
  cat(paste("Merging", length(objects_to_merge), "objects...\n"))
  merged_seurat <- merge(x = objects_to_merge[[1]], 
                         y = objects_to_merge[-1], 
                         add.cell.ids = names(objects_to_merge))
}

# Save output data
saveRDS(merged_seurat, file = output_path)
