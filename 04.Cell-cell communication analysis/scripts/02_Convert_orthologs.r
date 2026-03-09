#!/usr/bin/env Rscript
# Usage: Rscript 02_Convert_orthologs.R <merged.rds> <human_paralog.tsv> <species_paralog.tsv> <output_converted.rds>

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Error: Exactly 4 arguments required: <merged.rds> <human_paralog.tsv> <species_paralog.tsv> <output.rds>")
}

merged_rds_file <- args[1]
human_paralog_file <- args[2]
species_paralog_file <- args[3]
output_converted_rds_file <- args[4]

convert_genes_in_seurat <- function(seurat_rds_path, human_paralog_path, species_paralog_path, output_rds_path) {
  
  # Creating species-to-human gene map
  human_paralogs <- read.delim(human_paralog_path)
  species_paralogs <- read.delim(species_paralog_path)
  
  ortho_to_human_map <- human_paralogs %>% distinct(Orthogroup, .keep_all = TRUE) %>% select(Orthogroup, human_gene = Gene)
  species_to_ortho_map <- species_paralogs %>% distinct(Gene, .keep_all = TRUE) %>% select(species_gene = Gene, Orthogroup)
  species_to_human_df <- left_join(species_to_ortho_map, ortho_to_human_map, by = "Orthogroup") %>% filter(!is.na(human_gene))
  
  gene_map_vector <- setNames(species_to_human_df$human_gene, species_to_human_df$species_gene)
  
  # Loading Seurat object and extracting counts matrix
  seurat_obj <- readRDS(seurat_rds_path)
  counts_matrix <- seurat_obj[["RNA"]]$counts
  
  # Converting gene names and aggregating expression values
  original_genes <- rownames(counts_matrix)
  new_gene_names <- ifelse(original_genes %in% names(gene_map_vector), gene_map_vector[original_genes], original_genes)
  
  valid_indices <- !is.na(new_gene_names)
  counts_matrix <- counts_matrix[valid_indices, ]
  new_gene_names <- new_gene_names[valid_indices]
  original_genes <- original_genes[valid_indices]
  
  unique_new_names <- sort(unique(new_gene_names))
  
  mapping_matrix <- sparseMatrix(
    i = match(new_gene_names, unique_new_names),
    j = 1:length(new_gene_names),
    x = 1,
    dims = c(length(unique_new_names), length(new_gene_names)),
    dimnames = list(unique_new_names, original_genes)
  )
  
  aggregated_matrix <- mapping_matrix %*% counts_matrix
  
  # Creating and normalizing the new Seurat object
  new_seurat_obj <- CreateSeuratObject(counts = aggregated_matrix, meta.data = seurat_obj@meta.data)
  new_seurat_obj <- NormalizeData(new_seurat_obj)

  # Saving the output Seurat object
  saveRDS(new_seurat_obj, file = output_rds_path)
}

convert_genes_in_seurat(merged_rds_file, human_paralog_file, species_paralog_file, output_converted_rds_file)