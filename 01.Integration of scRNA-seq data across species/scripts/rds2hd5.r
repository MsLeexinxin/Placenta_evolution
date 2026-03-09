library(dplyr)
library(hdf5r)
library(Seurat)
library(SeuratDisk)
args<-commandArgs(TRUE)
infile<-args[1]
outname<-args[2]
data_fil<-readRDS(infile)
head(data_fil@meta.data)
levels(data_fil)
data_fil<-UpdateSeuratObject(data_fil)
raw_counts=GetAssayData(object =data_fil[['RNA']], slot = "counts")
data_fil[["RAW"]] <-CreateAssayObject(counts = raw_counts)
DefaultAssay(object = data_fil) <- "RAW"
SaveH5Seurat(data_fil, filename = paste(outname,".non_anno.h5Seurat",sep=""))
Convert(paste(outname,".non_anno.h5Seurat",sep=""), dest = "h5ad")

