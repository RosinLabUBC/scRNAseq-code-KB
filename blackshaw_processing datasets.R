# 4 data sets from 2 papers from E15 and E15.5 

library(Seurat)
library(dplyr)
library(patchwork)
library(clustree)


# uploading, quality control, SCTransform initial objects -----------------

# 2025 E15.5 Rep 1
file_path <- "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/blackshaw matrices/GSM8685159_TK67_filtered_feature_bc_matrix.h5"
Rep1_2025_data <- Read10X_h5(filename = file_path)
Rep1_2025 <- CreateSeuratObject(counts = Rep1_2025_data, project = "2025_E15.5_Rep1", min.cells = 3, min.features = 200)
Rep1_2025[["percent.mt"]] <- PercentageFeatureSet(Rep1_2025, pattern = "^mt-")
Rep1_2025 <- subset(Rep1_2025, subset = nFeature_RNA > 1000 & percent.mt < 5)

# 2025 E15.5 Rep 2
file_path_2 <- "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/blackshaw matrices/GSM8685160_TK68_filtered_feature_bc_matrix.h5"
Rep2_2025_data <- Read10X_h5(filename = file_path_2)
Rep2_2025 <- CreateSeuratObject(counts = Rep2_2025_data, project = "2025_E15.5_Rep2", min.cells = 3, min.features = 200)
Rep2_2025[["percent.mt"]] <- PercentageFeatureSet(Rep2_2025, pattern = "^mt-")
Rep2_2025 <- subset(Rep2_2025, subset = nFeature_RNA > 1000 & percent.mt < 5)


# 2020 E15 Rep 1
dir_path <- "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/blackshaw matrices/2020_e15_rep1"
Rep1_2020_data <- Read10X(data.dir = dir_path)
Rep1_2020 <- CreateSeuratObject(counts = Rep1_2020_data, project = "2020_E15_Rep1", min.cells = 3, min.features = 200)
Rep1_2020[["percent.mt"]] <- PercentageFeatureSet(Rep1_2020, pattern = "^mt-")
Rep1_2020 <- subset(Rep1_2020, subset = nFeature_RNA > 1000 & percent.mt < 5)


# 2020 E15 Rep 2
dir_path_rep2 <- "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/blackshaw matrices/2020_e15_rep2"
Rep2_2020_data <- Read10X(data.dir = dir_path_rep2)
Rep2_2020 <- CreateSeuratObject(counts = Rep2_2020_data, project = "2020_E15_Rep2", min.cells = 3, min.features = 200)
Rep2_2020[["percent.mt"]] <- PercentageFeatureSet(Rep2_2020, pattern = "^mt-")
Rep2_2020 <- subset(Rep2_2020, subset = nFeature_RNA > 1000 & percent.mt < 5)

# Set to 2 GiB
options(future.globals.maxSize = 2 * 1024^3)  

# SCTransform
Rep1_2025 <- SCTransform(Rep1_2025, assay = "RNA")
Rep2_2025 <- SCTransform(Rep2_2025, assay = "RNA")
Rep1_2020 <- SCTransform(Rep1_2020, assay = "RNA")
# Rep2_2020 <- SCTransform(Rep2_2020, assay = "RNA") -- not enough cells, going forward without this dataset (only 2 cells present after quality control)

# Set default assay
DefaultAssay(Rep1_2025) <- "SCT"
Rep1_2025 <- RunPCA(Rep1_2025)

# number of cells
n_cells_Rep1_2025 <- ncol(Rep1_2025)
print(paste("2025 E15.5 Rep 1 cells:", n_cells_Rep1_2025))
# 5281 cells

# Save object
saveRDS(Rep1_2025, "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/saved objects/Rep1_2025.rds")


# Set default assay
DefaultAssay(Rep2_2025) <- "SCT"
Rep2_2025 <- RunPCA(Rep2_2025)

# number of cells
n_cells_Rep2_2025 <- ncol(Rep2_2025)
print(paste("2025 E15.5 Rep 2 cells:", n_cells_Rep2_2025))
# 6112 cells

# Save object
saveRDS(Rep2_2025, "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/saved objects/Rep2_2025.rds")


# Set default assay
DefaultAssay(Rep1_2020) <- "SCT"
Rep1_2020 <- RunPCA(Rep1_2020)

# number of cells
n_cells_Rep1_2020 <- ncol(Rep1_2020)
print(paste("2020 E15 Rep 1 cells:", n_cells_Rep1_2020))
# 292 cells

# Save object
saveRDS(Rep1_2020, "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/saved objects/Rep1_2020.rds")

# merging objects  -------------------------------------------------

# load in objects 
Rep1_2025 <- readRDS('/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/saved objects/Rep1_2025.rds')
Rep2_2025 <- readRDS('/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/saved objects/Rep2_2025.rds')
Rep1_2020 <- readRDS('/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/saved objects/Rep1_2020.rds')

# Merge COLD_NSC and CON_NSC
merged_object <- merge(
  x = Rep1_2025,                          # reference object
  y = c(Rep2_2025, Rep1_2020),            # object to merge
  add.cell.ids = c("2025_R1", "2025_R2", "2020_R1"),     # prefixes to distinguish cells
  project = "Merged_NSC")


saveRDS(merged_object, "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/saved objects/merged_object.rds")

# checked umap on merged object before SCTransforming and the cells were distributed nicely among eachother, don't need to integrate 

# analysis without integration --------------------------------------------

# load in merged_object
merged_object <- readRDS('/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/saved objects/merged_object.rds')

# Set to 2 GiB
options(future.globals.maxSize = 2 * 1024^3)  

# SCTransform merged object
merged_object <- SCTransform(object = merged_object, assay = "RNA")

DefaultAssay(merged_object) <- "SCT"
merged_object <- RunPCA(merged_object)
ElbowPlot(merged_object, ndims = 40, reduction = "pca")		

saveRDS(merged_object, "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/saved objects/merged_object_SCTransformed.rds")

# clustering merged object that is SCTransformed --------------------------
# merged object contains 3 blackshaw CON data sets, they have been SCTransformed again

# load in merged_object
merged_object <- readRDS('/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/saved objects/merged_object_SCTransformed.rds')

# Confirm number of cells
ncol(merged_object)
# 11685 cells

# pull out Vim+ cells from merged object ----------------------------------

# Subset cells expressing vimentin (Vim > 0)
merged_object_vim_pos <- subset(merged_object, subset = Vim > 0)

# Confirm number of cells
ncol(merged_object_vim_pos)
# 5571 cells

DefaultAssay(merged_object_vim_pos) <- "SCT"
merged_object_vim_pos <- RunPCA(merged_object_vim_pos)
ElbowPlot(merged_object_vim_pos, ndims = 40, reduction = "pca")		

saveRDS(merged_object_vim_pos, "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/saved objects/merged_object_vim_pos.rds")

# load in Vim+ subsetted cells --------------------------------------------
merged_object_vim_pos <- readRDS('/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/saved objects/merged_object_vim_pos.rds')

# Confirm number of cells
ncol(merged_object_vim_pos)
# 5571 cells

# using 35 dims, 0.05 res to cluster (found on separate workflow "testing parameters to cluster on Vim+.R")
# 0.05 res looked the best (33 dims, 0.05 res was similar but going forward with 35 dims, 0.05 res to make it more similar to previous workflow)
merged_object_vim_pos <- FindNeighbors(merged_object_vim_pos, dims = 1:35)	
merged_object_vim_pos <- FindClusters(merged_object_vim_pos, resolution = 0.05)		
merged_object_vim_pos <- RunUMAP(merged_object_vim_pos, dims = 1:35)
DimPlot(merged_object_vim_pos, reduction = "umap")
DimPlot(merged_object_vim_pos, reduction = "umap", group.by = c("orig.ident", "seurat_clusters"), split.by = "orig.ident")

# subset out cluster 1 and Ascl1+ cells ----------------------------------------------------
# Using WhichCells() + union() ensures each cell appears only once in the final subset and because you are subsetting from the same object 

# Step 1: Subset cells from cluster 1
cluster1_cells <- WhichCells(merged_object_vim_pos, idents = "1")

# Step 2: Subset cells with Ascl1 > 0
ascl1_pos_cells <- WhichCells(merged_object_vim_pos, expression = Ascl1 > 0)

# Step 3: Combine the two sets of cell names (union, not intersection)
combined_cells <- union(cluster1_cells, ascl1_pos_cells)

# Step 4: Subset the original object using these cell names
merged_object_vim_pos_subset <- subset(merged_object_vim_pos, cells = combined_cells)

# Step 5: Check the original sample identities (should still be CON and COLD)
table(merged_object_vim_pos_subset$orig.ident)
#    2020_E15_Rep1   2025_E15.5_Rep1   2025_E15.5_Rep2 
#    56              1151              1340     (total = 2547)

DefaultAssay(merged_object_vim_pos_subset) <- "SCT"
merged_object_vim_pos_subset <- RunPCA(merged_object_vim_pos_subset)
ElbowPlot(merged_object_vim_pos_subset, ndims = 40, reduction = "pca")		

saveRDS(merged_object_vim_pos_subset, "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/saved objects/merged_object_vim_pos_subset.rds")


# load in Vim+ cells that have been subsetted (cluster 1 or Ascl1+) -----------------------------------------

merged_object_vim_pos_subset <- readRDS('/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/saved objects/merged_object_vim_pos_subset.rds')

# testing parameters on cluster 1 + Ascl1+ subset on separate workflow 

# using 33 dims, 0.2 res -----------------------------
merged_object_vim_pos_subset <- FindNeighbors(merged_object_vim_pos_subset, dims = 1:33)	
merged_object_vim_pos_subset <- FindClusters(merged_object_vim_pos_subset, resolution = 0.2)		
merged_object_vim_pos_subset <- RunUMAP(merged_object_vim_pos_subset, dims = 1:33)
DimPlot(merged_object_vim_pos_subset, reduction = "umap")
DimPlot(merged_object_vim_pos_subset, reduction = "umap", group.by = c("orig.ident", "seurat_clusters"), split.by = "orig.ident")

# quality check  ----------------------------------------------------------

# Total number of cells 
ncol(merged_object_vim_pos_subset)
# 2547

# Number of cells in each cluster 
table(Idents(merged_object_vim_pos_subset))
#    0    1    2    3    4    5 
#  901  553  405  359   244   85 

# Number of CON vs. COLD cells in total
table(merged_object_vim_pos_subset$orig.ident)
#     2020_E15_Rep1   2025_E15.5_Rep1   2025_E15.5_Rep2 
#      56              1151              1340 

# Number of CON and COLD cells in each cluster
table(merged_object_vim_pos_subset$seurat_clusters, merged_object_vim_pos_subset$orig.ident)
#          2020_E15_Rep1  2025_E15.5_Rep1  2025_E15.5_Rep2
# 0            14             416             471
# 1             2             240             311
# 2            38             177             190
# 3             0             170             189
# 4             2             110             132
# 5             0              38              47

# removing male and female double positive and non-expressing cells -------
# going to just use "Xist" for females and "Ddx3y", "Uty", "Eif2s3y" for males

# Xist – X-inactive specific transcript (a long non-coding RNA that inactivates one of the X chromosomes in females)

# Ddx3y – DEAD-box helicase 3 Y-linked
# Uty – Ubiquitously transcribed Y-linked
# Eif2s3y – Eukaryotic translation initiation factor 2 subunit 3 Y-linked

# define markers and expression threshold
female_genes <- c("Xist")
male_genes <- c("Ddx3y", "Uty", "Eif2s3y")
expr_threshold <- 0.1

# Extract SCTransformed data matrix
expr <- GetAssayData(merged_object_vim_pos_subset, assay = "SCT", layer = "data")

# Identify cells expressing male or female markers
female_expr <- colSums(expr[female_genes, , drop = FALSE] > expr_threshold) > 0
male_expr   <- colSums(expr[male_genes, , drop = FALSE] > expr_threshold) > 0

# Assign sex labels: only keep non-overlapping classifications
sex_labels <- ifelse(female_expr & !male_expr, "Female",
                     ifelse(male_expr & !female_expr, "Male", NA))

# Store in metadata
merged_object_vim_pos_subset$sex <- sex_labels

# Directly use the metadata to filter valid cells
valid_cells <- rownames(merged_object_vim_pos_subset@meta.data)[!is.na(merged_object_vim_pos_subset@meta.data$sex)]

# Subset by cell names
clean_combined_object <- subset(merged_object_vim_pos_subset, cells = valid_cells)

# Optional check
table(clean_combined_object$sex)
# Female   Male 
# 593     1270    (Total 1863)

saveRDS(clean_combined_object, file = "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/saved objects/clean_combined_object_blackshaw.rds")


# identifying Cycling NSPC cells ------------------------------------------

# cell cycle analysis identified cluster 3 as cycling (done on 'cell cycle workflow on blackshaw cells.R' workflow)

# Cycling NSPCs cell markers
plot <- FeaturePlot(clean_combined_object, features = c("Birc5",
                                                        "Ube2c",
                                                        "Ccnb1",
                                                        "Cdca8",
                                                        "Tpx2"))
ggsave(filename = '/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/analysis on Vim+ subsetted cells/Cycling NSPCs.png', plot, width = 15, height = 15)

# cluster 3 is cycling and has cycling NSPC markers 

# naming cluster 3  -------------------------------------------------------
clean_combined_object$cluster_names <- ifelse(
  as.character(clean_combined_object$seurat_clusters) == "3", 
  "Cycling NSPCs", 
  as.character(clean_combined_object$seurat_clusters))

# colour
unique_names <- unique(clean_combined_object$cluster_names)
cluster_colors1 <- scales::hue_pal()(length(unique_names))
names(cluster_colors1) <- unique_names

# Overwrite the "Cycling NSPCs" slot with your specific green color
cluster_colors1["Cycling NSPCs"] <- "#00BA38"

# cluster 0 is too close to cycling nspcs
cluster_colors1["0"] <- "#FF7F50"

# store in metadata
clean_combined_object$treatment <- "CON"
clean_combined_object$replicate_group <- clean_combined_object$orig.ident

# Plot UMAP using your custom names and colors
DimPlot(
  clean_combined_object, 
  reduction = "umap", 
  group.by = "cluster_names", 
  cols = cluster_colors1
) + ggtitle("Merged Blackshaw NSC Dataset - Replicate Review")

# Print tables to verify cell distributions across your different replicate groups
table(clean_combined_object$replicate_group)
#   2020_E15_Rep1   2025_E15.5_Rep1   2025_E15.5_Rep2 
#       28               829             1006 

table(clean_combined_object$cluster_names)
#        0             1             2             4             5      Cycling NSPCs 
#      640           415           297           171            66           274 


saveRDS(clean_combined_object, file = "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/saved objects/clean_combined_object_blackshaw_classified_cycling_nspcs.rds")


# load in clean combined object -------------------------------------------

clean_combined_object <- readRDS('/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/saved objects/clean_combined_object_blackshaw_classified_cycling_nspcs.rds')

# This creates a 2-way table showing the male/female split for ALL clusters
table(clean_combined_object$cluster_names, clean_combined_object$sex)

#                 Female Male
# 0                191  449
# 1                113  302
# 2                123  174
# 4                 53  118
# 5                 15   51
# Cycling NSPCs     98  176

