library(Seurat)
library(muscat)
library(SingleCellExperiment)
library(dplyr)
library(tibble)
library(data.table)

args <- commandArgs(TRUE)
infile  <- args[1]   # input Seurat RDS
outname <- args[2]   # output prefix

seu <- readRDS(infile)
DefaultAssay(seu) <- "RNA"
sce <- as.SingleCellExperiment(seu)
sce <- prepSCE(
  sce,
  kid  = "seurat_clusters",
  gid  = "cell_type",
  sid  = "batch",
  drop = TRUE
)
pb <- aggregateData(
  sce,
  assay = "counts",
  fun   = "sum",
  by    = c("group_id")
)

pb_counts <- pb@assays@data@listData[[1]]
pb_counts_df <- pb_counts %>%
  as.data.frame() %>%
  rownames_to_column(var = "Gene")

lib_size <- colSums(pb_counts)
pb_cpm <- t( t(pb_counts) / lib_size * 1e6 )

pb_cpm_df <- pb_cpm %>%
  as.data.frame() %>%
  rownames_to_column(var = "Gene")

fwrite(
  pb_cpm_df,
  file = paste0(outname, ".pseudobulk.CPM.txt"),
  sep  = "\t"
)
