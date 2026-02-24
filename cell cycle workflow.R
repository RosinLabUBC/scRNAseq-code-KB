# https://satijalab.org/seurat/archive/v4.3/cell_cycle_vignette , https://github.com/satijalab/seurat/blob/HEAD/vignettes/cell_cycle_vignette.Rmd 

# Running cell cycle analysis on all cells (not cleaned), cleaned merged cells, clean female cells, clean male cells

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
clean_combined_object <- readRDS('/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/clean_combined_object_CON&COLD_strat2.rds')

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

# on Female cells -------------------------------------
female_obj <- readRDS('/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/clean_female_obj_CON&COLD_strat2.rds')

# Run cell cycle scoring using mouse gene symbols
female_obj <- CellCycleScoring(
  female_obj,
  s.features = mmus_s,
  g2m.features = mmus_g2m,
  set.ident = TRUE)

# View scores and predicted phase
head(female_obj[[]])

# DimPlot by predicted cell cycle phase (from CellCycleScoring)
DimPlot(female_obj, group.by = "Phase", reduction = "umap") + ggtitle("Cell Cycle Phase before Regression_Clean Female Only")

# on Male cells -------------------------------------
male_obj <- readRDS('/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/clean_male_obj_CON&COLD_strat2.rds')

# Run cell cycle scoring using mouse gene symbols
male_obj <- CellCycleScoring(
  male_obj,
  s.features = mmus_s,
  g2m.features = mmus_g2m,
  set.ident = TRUE)

# View scores and predicted phase
head(male_obj[[]])

# DimPlot by predicted cell cycle phase (from CellCycleScoring)
DimPlot(male_obj, group.by = "Phase", reduction = "umap") + ggtitle("Cell Cycle Phase before Regression_Clean Male Only")


# on CON cells -------------------------------------
con_obj <- readRDS('/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/clean_con_obj_M&F_strat2.rds')

# Run cell cycle scoring using mouse gene symbols
con_obj <- CellCycleScoring(
  con_obj,
  s.features = mmus_s,
  g2m.features = mmus_g2m,
  set.ident = TRUE)

# View scores and predicted phase
head(con_obj[[]])

# DimPlot by predicted cell cycle phase (from CellCycleScoring)
DimPlot(con_obj, group.by = "Phase", reduction = "umap") + ggtitle("Cell Cycle Phase before Regression_Clean CON Only")


# on COLD cells -------------------------------------
cold_obj <- readRDS('/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/clean_cold_obj_M&F_strat2.rds')

# Run cell cycle scoring using mouse gene symbols
cold_obj <- CellCycleScoring(
  cold_obj,
  s.features = mmus_s,
  g2m.features = mmus_g2m,
  set.ident = TRUE)

# View scores and predicted phase
head(cold_obj[[]])

# DimPlot by predicted cell cycle phase (from CellCycleScoring)
DimPlot(cold_obj, group.by = "Phase", reduction = "umap") + ggtitle("Cell Cycle Phase before Regression_Clean COLD Only")
