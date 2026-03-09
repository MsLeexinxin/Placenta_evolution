library(Seurat)
library(dplyr)
library(harmony)
args<-commandArgs(TRUE)
control_file<-args[1]
case_file<-args[2]
outname<-args[3]

control <- readRDS(control_file)
case <-readRDS(case_file)

control$batch <- "control"
case$batch <- "case"

if (class(control$RNA) == 'Assay5'){control$RNA = as(control$RNA, 'Assay')}
if (class(case$RNA) == 'Assay5'){case$RNA = as(case$RNA, 'Assay')}

combined <- merge(control, y = list(case), add.cell.ids = c("CT","CASE"))
combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)
combined <- ScaleData(combined)
combined <- RunPCA(combined, features = VariableFeatures(object = combined))
combined <- RunHarmony(combined, group.by.vars = "batch")


control_others <- FindMarkers(combined,group.by = "batch", ident.1="control",min.pct = 0.1,testuse = "wilcox",logfc.threshold = 0.5)
write.csv(control_others,paste(outname,".control.DEG.harmony.csv",sep=""))

case_others <- FindMarkers(combined,group.by = "batch", ident.1="case",min.pct = 0.1,testuse = "wilcox",logfc.threshold = 0.5)
write.csv(case_others,paste(outname,".case_others.DEG.harmony.csv",sep=""))


saveRDS(combined,paste(outname,".merge.harmony.rds",sep=""))

