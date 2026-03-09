library(Seurat)
library(dplyr)
library(harmony)
args<-commandArgs(TRUE)
h_file<-args[1]
ma_file<-args[2]
m_file<-args[3]
t_file<-args[4]
d_file<-args[5]
s_file<-args[6]
g_file<-args[7]
p_file<-args[8]
homolog_file<-args[9]
outname<-args[10]

human <- readRDS(h_file)
macaca<-readRDS(ma_file)
mouse<-readRDS(m_file)
treeshrew<-readRDS(t_file)
dog<-readRDS(d_file)
sheep<-readRDS(s_file)
goat<-readRDS(g_file)
pig<-readRDS(p_file)
homolog_genes <- read.csv(homolog_file, header = TRUE)

human_genes <- intersect(homolog_genes$Human,rownames(human))
macaca_genes <- intersect(homolog_genes$Macaca,rownames(macaca))
mouse_genes <- intersect(homolog_genes$Mouse,rownames(mouse))
treeshrew_genes <- intersect(homolog_genes$Treeshrew,rownames(treeshrew))
dog_genes <- intersect(homolog_genes$Dog,rownames(dog))
sheep_genes <- intersect(homolog_genes$Sheep,rownames(sheep))
goat_genes <- intersect(homolog_genes$Goat,rownames(goat))
pig_genes <- intersect(homolog_genes$Pig,rownames(pig))

human_org<-homolog_genes[match(human_genes,homolog_genes$Human),1]
macaca_org<-homolog_genes[match(macaca_genes,homolog_genes$Macaca),1]
mouse_org<-homolog_genes[match(mouse_genes,homolog_genes$Mouse),1]
treeshrew_org<-homolog_genes[match(treeshrew_genes,homolog_genes$Treeshrew),1]
dog_org<-homolog_genes[match(dog_genes,homolog_genes$Dog),1]
sheep_org<-homolog_genes[match(sheep_genes,homolog_genes$Sheep),1]
goat_org<-homolog_genes[match(goat_genes,homolog_genes$Goat),1]
pig_org<-homolog_genes[match(pig_genes,homolog_genes$Pig),1]

gene_lists <- list(human_org,macaca_org,mouse_org,treeshrew_org,dog_org,sheep_org,goat_org,pig_org)

common_org <- Reduce(intersect, gene_lists)

common_genes <- homolog_genes[homolog_genes$Orthogroup %in% common_org, ]

human_sub <- subset(human, features = common_genes$Human)
macaca_sub <- subset(macaca, features = common_genes$Macaca)
mouse_sub <- subset(mouse, features = common_genes$Mouse)
treeshrew_sub <- subset(treeshrew, features = common_genes$Treeshrew)
dog_sub <- subset(dog, features = common_genes$Dog)
sheep_sub <- subset(sheep, features = common_genes$Sheep)
goat_sub <- subset(goat, features = common_genes$Goat)
pig_sub <- subset(pig, features = common_genes$Pig)

#human@assays$RNA@counts <- human@assays$RNA@counts[common_genes$Human,]
#macaca@assays$RNA@counts <- macaca@assays$RNA@counts[common_genes$Human,]
#mouse@assays$RNA@counts <- mouse@assays$RNA@counts[common_genes$Human,]
#treeshrew@assays$RNA@counts <- treeshrew@assays$RNA@counts[common_genes$Human,]
#dog@assays$RNA@counts <- dog@assays$RNA@counts[common_genes$Human,]
#sheep@assays$RNA@counts  <- sheep@assays$RNA@counts[common_genes$Human,]
#goat@assays$RNA@counts <-goat@assays$RNA@counts[common_genes$Human,]
#pig@assays$RNA@counts<- pig@assays$RNA@counts[common_genes$Human,]

if (class(macaca_sub$RNA) == 'Assay5'){macaca_sub$RNA = as(macaca_sub$RNA, 'Assay')}

rownames(human_sub@assays$RNA@counts) <- common_genes$Human
rownames(macaca_sub@assays$RNA@counts) <- common_genes$Human
rownames(mouse_sub@assays$RNA@counts) <- common_genes$Human
rownames(treeshrew_sub@assays$RNA@counts) <- common_genes$Human
rownames(dog_sub@assays$RNA@counts) <- common_genes$Human
rownames(sheep_sub@assays$RNA@counts) <- common_genes$Human
rownames(goat_sub@assays$RNA@counts) <- common_genes$Human
rownames(pig_sub@assays$RNA@counts) <- common_genes$Human

rownames(human_sub@assays$RNA@data) <- common_genes$Human
rownames(macaca_sub@assays$RNA@data) <- common_genes$Human
rownames(mouse_sub@assays$RNA@data) <- common_genes$Human
rownames(treeshrew_sub@assays$RNA@data) <- common_genes$Human
rownames(dog_sub@assays$RNA@data) <- common_genes$Human
rownames(sheep_sub@assays$RNA@data) <- common_genes$Human
rownames(goat_sub@assays$RNA@data) <- common_genes$Human
rownames(pig_sub@assays$RNA@data) <- common_genes$Human

human_sub$batch <- "human"
macaca_sub$batch <- "macaca"
mouse_sub$batch <- "mouse"
treeshrew_sub$batch <- "treeshrew"
dog_sub$batch <- "dog"
sheep_sub$batch <- "sheep"
goat_sub$batch <- "goat"
pig_sub$batch <- "pig"

combined <- merge(human_sub, y = list(macaca_sub,mouse_sub,treeshrew_sub,dog_sub,sheep_sub,goat_sub,pig_sub), add.cell.ids = c("HO","MA", "MU", "TU","DO","OV","GO","SU"))
combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)
combined <- ScaleData(combined)
combined <- RunPCA(combined, features = VariableFeatures(object = combined))
combined <- RunHarmony(combined, group.by.vars = "batch")

combined@meta.data <- combined@meta.data %>% mutate(group = case_when(batch %in% c("human", "macaca", "mouse") ~ "Hemo",batch %in% c("dog", "treeshrew") ~ "Endo",batch %in% c("sheep", "goat","pig") ~ "Epi",TRUE ~ NA))

hemo_others <- FindMarkers(combined,group.by = "group", ident.1="Hemo",min.pct = 0.1,testuse = "wilcox",logfc.threshold = 0.5)
write.csv(hemo_others,paste(outname,".hemo_others.DEG.harmony.csv",sep=""))

endo_others <- FindMarkers(combined,group.by = "group", ident.1="Endo",min.pct = 0.1,testuse = "wilcox",logfc.threshold = 0.5)
write.csv(endo_others,paste(outname,".endo_others.DEG.harmony.csv",sep=""))

epi_others<-FindMarkers(combined,group.by = "group", ident.1="Epi",min.pct = 0.1,testuse = "wilcox",logfc.threshold = 0.5)
write.csv(epi_others,paste(outname,".epi_others.DEG.harmony.csv",sep=""))

saveRDS(combined,paste(outname,".merge.harmony.rds",sep=""))

