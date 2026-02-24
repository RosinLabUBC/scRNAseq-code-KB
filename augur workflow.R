# https://github.com/neurorestore/Augur 
# for DEG Analysis and Augur Analysis, switch to RNA Assay and remember to Normalize, FindVariableFeatures, ScaleData
# Done in "MERGED_CON&COLD.R" and "Vim expressing cells.R" to get clean_female_merged_object
# removed cells with both M and F markers and cells with no sex markers from the same object, (didn't split by sex)
# splitting by sex after removing unwanted cells 

library(SeuratExtend)
library(Seurat)
library(dplyr)
library(patchwork)
library(clustree)
library(ggplot2)
library(Augur)


# load in female and male objects -----------------------------------------

female_obj <- readRDS('/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/clean_female_obj_CON&COLD_strat2.rds')
male_obj <- readRDS('/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/clean_male_obj_CON&COLD_strat2.rds')

# Set default assay to RNA
DefaultAssay(female_obj) <- "RNA"
DefaultAssay(male_obj) <- "RNA"

# Normalize and scale RNA assay
female_obj <- NormalizeData(female_obj)
female_obj <- FindVariableFeatures(female_obj)
female_obj <- ScaleData(female_obj, features = row.names(female_obj))

male_obj <- NormalizeData(male_obj)
male_obj <- FindVariableFeatures(male_obj)
male_obj <- ScaleData(male_obj, features = row.names(male_obj))

# JoinLayers to properly link for analysis 
female_obj <- JoinLayers(female_obj, layers = "data")
male_obj <- JoinLayers(male_obj, layers = "data")


# Female only, default settings -------------------------------------------

# Set the seed — ensures reproducibility of Augur's random subsampling
set.seed(102)

augur_result_female <- calculate_auc(
  female_obj,
  label_col = "treatment",               # COLD vs CON
  cell_type_col = "cluster_names",     # cluster labels
  n_threads = 8) # n_threads = 8 uses 8 of your 14 available cores → good balance between speed and system responsiveness

# Save AUC scores
write.csv(
  augur_result_female$AUC,
  "/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Vim+ and subsetted cells/Augur/female_augur_auc_scores(labelled).csv",
  row.names = FALSE)

# Plot
ggplot(augur_result_female$AUC, aes(x = reorder(cell_type, auc), y = auc)) +
  geom_col(fill = "palevioletred2") +
  coord_flip() +
  labs(x = "Cluster", y = "AUC (Perturbation Score)", title = "Augur Prioritization (Female Cells)") +
  theme_minimal()

# Plot (formatted for thesis figure)
ggplot(augur_result_female$AUC, aes(x = reorder(cell_type, auc), y = auc)) +
  geom_col(fill = "palevioletred2", width = 0.5) +  # smaller width → skinnier bars
  coord_flip() +
  labs(
    x = "Cluster",
    y = "AUC (Perturbation Score)",
    title = "Augur Prioritization (Female Cells)") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14))

# plot umap (rank)
plot_umap(
  augur = augur_result_female,
  sc = female_obj,
  mode = "rank",              # Show ranked AUCs (relative prioritization across cell types)
  reduction = "umap",
  palette = c("blue", "gray", "red"),  # blue = low rank (low perturbation), red = high
  cell_type_col = "cluster_names")

# plot umap (raw auc)
plot_umap(
  augur = augur_result_female,
  sc = female_obj,
  mode = "default",              # Show raw AUC values
  reduction = "umap",
  palette = c("blue", "gray", "red"),  # blue = low rank (low perturbation), red = high
  cell_type_col = "cluster_names")


# Male only, default settings -------------------------------------------

set.seed(103)

augur_result_male <- calculate_auc(
  male_obj,
  label_col = "treatment",
  cell_type_col = "cluster_names",
  n_threads = 8) # n_threads = 8 uses 8 of your 14 available cores → good balance between speed and system responsiveness 

# Save AUC scores
write.csv(
  augur_result_male$AUC,
  "/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Vim+ and subsetted cells/Augur/male_augur_auc_scores(labelled).csv",
  row.names = FALSE)

# Plot
ggplot(augur_result_male$AUC, aes(x = reorder(cell_type, auc), y = auc)) +
  geom_col(fill = "steelblue3") +
  coord_flip() +
  labs(x = "Cluster", y = "AUC (Perturbation Score)", title = "Augur Prioritization (Male Cells)") +
  theme_minimal()

# Plot (formatted for thesis figure)
ggplot(augur_result_male$AUC, aes(x = reorder(cell_type, auc), y = auc)) +
  geom_col(fill = "steelblue3", width = 0.5) +  # smaller width → skinnier bars
  coord_flip() +
  labs(
    x = "Cluster",
    y = "AUC (Perturbation Score)",
    title = "Augur Prioritization (Male Cells)") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14))

# plot umap (rank)
plot_umap(
  augur = augur_result_male,
  sc = male_obj,
  mode = "rank",              # Show ranked AUCs (relative prioritization across cell types)
  reduction = "umap",
  palette = c("blue", "gray", "red"),  # blue = low rank (low perturbation), red = high
  cell_type_col = "cluster_names")

# plot umap (raw auc)
plot_umap(
  augur = augur_result_male,
  sc = male_obj,
  mode = "default",              # Show raw AUC values
  reduction = "umap",
  palette = c("blue", "gray", "red"),  # blue = low rank (low perturbation), red = high
  cell_type_col = "cluster_names")


