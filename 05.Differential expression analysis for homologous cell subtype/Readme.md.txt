Cross-Species Placental Single-Cell Analysis Pipeline for Differential Gene Expression Across Placental Types

This README outlines a pipeline for analyzing single-cell RNA-seq data across eight mammalian species to identify differentially expressed genes associated with different placental types. The pipeline integrates data from species with different placental morphologies and performs comparative differential expression analysis.

Overview

This pipeline integrates single-cell expression data from eight mammalian species, groups them by placental type, and performs differential gene expression analysis to identify genes associated with hemochorial, endotheliochorial, and epitheliochorial placentation.

Species and Data Files

Input RDS Files

The pipeline requires the following single-cell RDS files, each representing the same cell type across different species:

File Name Species Code

H1_H2.anno.rds Human H

MA1_MA2_new.anno.rds Macaque (Macaca) MA

M1_M2.anno.rds Mouse M

SQ1_SQ2.anno.rds Treeshrew SQ

D1_D2.anno.rds Dog D

S1_S2.anno.rds Sheep S

G1_G2.anno.rds Goat G

Z1_Z2.anno.rds Pig Z

Orthologous Genes File

H_MA_M_SQ_D_S_G_Z.org.csv - A CSV file containing orthologous gene groups across all eight species.

File Format:

Orthogroup,Human,Macaca,Mouse,Treeshrew,Dog,Sheep,Goat,Pig
PRDX1,PRDX1,PRDX1,Prdx1,PRDX1,PRDX1,PRDX1,PRDX1,PRDX1
RPL35,RPL35,RPL35,Rpl35,RPL35,RPL35,RPL35,RPL35,RPL35
CERS4,CERS4,CERS4,Cers4,CERS4,CERS4,CERS4,CERS4,CERS4
SRA1,SRA1,SRA1,Sra1,SRA1,SRA1,SRA1,SRA1,SRA1
PARP14,PARP14,PARP14,Parp14,PARP14,PARP14,PARP14,PARP14,PARP14


Columns:
• Orthogroup: Gene family/orthogroup identifier

• Human: Human gene symbol

• Macaca: Macaque gene symbol

• Mouse: Mouse gene symbol

• Treeshrew: Treeshrew gene symbol

• Dog: Dog gene symbol

• Sheep: Sheep gene symbol

• Goat: Goat gene symbol

• Pig: Pig gene symbol

Placental Type Classification

Species are grouped into three placental types based on their placental morphology:

Placental Type Species Included Species

Hemochorial Hemo Human, Macaque, Mouse

Endotheliochorial Endo Treeshrew, Dog

Epitheliochorial Epi Sheep, Goat, Pig

Pipeline Execution

Main Analysis Script

Run the following command to perform the cross-species integration and differential expression analysis:
Rscript combine.cluster.R \
  H1_H2.anno.rds \
  MA1_MA2_new.anno.rds \
  M1_M2.anno.rds \
  SQ1_SQ2.anno.rds \
  D1_D2.anno.rds \
  S1_S2.anno.rds \
  G1_G2.anno.rds \
  Z1_Z2.anno.rds \
  H_MA_M_SQ_D_S_G_Z.org.csv \
  8species_Tro


Parameters:
1. Eight RDS files: Single-cell data files for each species
2. Orthologous genes file: CSV file with gene mappings across species
3. Output prefix: 8species_Tro (customizable output file prefix)

What the Script Does

The combine.cluster.R script performs the following steps:

1. Data Integration: Loads and integrates single-cell data from all eight species
2. Harmony Integration: Uses Harmony algorithm to integrate data and remove batch effects
3. Placental Type Grouping: Groups species into three placental types as described above
4. Differential Expression Analysis: Performs pairwise comparisons between placental types using the Wilcoxon rank-sum test (test.use = "wilcox") via the FindMarkers function
5. Output Generation: Produces the following output files

Output Files

The script generates four main output files with the specified prefix:

1. 8species_Tro.hemo_others.DEG.harmony.csv
   • Contains differentially expressed genes (DEGs) when comparing Hemochorial placenta against the other two types

   • Columns include gene identifiers, average log2 fold changes, p-values, adjusted p-values, and expression statistics

2. 8species_Tro.endo_others.DEG.harmony.csv
   • Contains DEGs when comparing Endotheliochorial placenta against the other two types

3. 8species_Tro.epi_others.DEG.harmony.csv
   • Contains DEGs when comparing Epitheliochorial placenta against the other two types

4. 8species_Tro.merge.harmony.rds
   • The integrated Seurat object containing all eight species' data after Harmony batch correction

   • Can be loaded in R for further visualization and analysis: seurat_obj <- readRDS("8species_Tro.merge.harmony.rds")

Output File Structure

The CSV output files contain the following columns (example from DESeq2/Wilcoxon test output):
• gene: Gene identifier

• avg_log2FC: Average log2 fold change

• pct.1: Percentage of cells expressing the gene in group 1

• pct.2: Percentage of cells expressing the gene in group 2

• p_val: Raw p-value

• p_val_adj: Adjusted p-value (Bonferroni or similar correction)

• cluster: Cell cluster information

• gene_symbol: Gene symbol (if available)

Disease Case-Control Analysis Extension

For disease case-control comparisons (e.g., comparing disease vs. control samples), a modified version of the script is available:
Rscript combine.cluster.disease.R control.rds case.rds outname


Parameters:
• control.rds: Control group single-cell data file

• case.rds: Disease/case group single-cell data file

• outname: Output file prefix

Output Files:
1. outname.control.DEG.harmony.csv: Upregulated genes in the control group compared to the case group
2. outname.case_others.DEG.harmony.csv: Upregulated genes in the case group compared to the control group
3. outname.merge.harmony.rds: Integrated Seurat object containing both control and case data after Harmony batch correction

Requirements

• R (version 4.0 or higher)

• Seurat package (for single-cell analysis)

• Harmony package (for data integration)

• Other dependencies: dplyr, ggplot2, etc.

Usage Notes

1. Ensure all input files are in the same directory or provide full paths
2. The orthologous genes file must contain exactly the eight species columns as shown
3. The RDS files should contain annotated Seurat objects with consistent cell type annotations
4. Adjust the output prefix (8species_Tro) as needed for different analyses
5. For disease analysis, ensure control and case samples are properly annotated in their respective RDS files

Interpretation

The output files can be used to:
1. Identify genes specifically upregulated in each placental type
2. Understand molecular differences between placental morphologies
3. Discover evolutionary adaptations in placental development
4. Generate hypotheses for functional validation studies
5. For disease studies, identify differentially expressed genes between disease and control conditions