library(Seurat)
library(dplyr)
library(patchwork)
library(clustree)
library(biomaRt)
library(gprofiler2)

# Load the Seurat human cell cycle genes
s.genes_human <- cc.genes.updated.2019$s.genes
g2m.genes_human <- cc.genes.updated.2019$g2m.genes

# Convert them to mouse using gorth()
mmus_s <- gorth(s.genes_human, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
mmus_g2m <- gorth(g2m.genes_human, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name

# on all cells, sex filtered -------------------------------------
clean_combined_object <- readRDS('/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/blackshaw + kb/clean_combined_object_blackshaw_and_kb.rds')

# Run cell cycle scoring using mouse gene symbols
clean_combined_object <- CellCycleScoring(
  clean_combined_object,
  s.features = mmus_s,
  g2m.features = mmus_g2m,
  set.ident = TRUE)

# View scores and predicted phase
head(clean_combined_object[[]])

# DimPlot by predicted cell cycle phase (from CellCycleScoring)
DimPlot(clean_combined_object, group.by = "Phase", reduction = "umap") + ggtitle("Cell Cycle Phase before Regression_Clean Combined")
