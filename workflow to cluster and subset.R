# On Vim+ cells only, pull out all of cluster 1 and Ascl1+ and cluster together from 35 dims, 0.1 res
# looks pretty similar to 35 dims, 0.05 res, but 0.1 res seems to have less red interspersed in cluster 1 (more cells in cluster 1)

# Subset out Vim+ cells
# Subset out cluster 1 (contains cycling/proliferative and RGCs) and Ascl1+ cells (which is a marker for NPCs)
# Subset out cells double positive for male and female sex markers and ones non-expressing
# Do analysis on combined object and male and female split from the combined object

library(Seurat)
library(dplyr)
library(patchwork)
library(clustree)

# COLD Cells (Male and Female) --------------------------------------------

expression_matrix <- ReadMtx(mtx = "/Users/kodeebao/Documents/StephaneMatrix/matrices/Apr2025_cold_rna_velocity/counts_unfiltered/cells_x_genes.cell.mtx", features = "/Users/kodeebao/Documents/StephaneMatrix/matrices/Apr2025_cold_rna_velocity/counts_unfiltered/cells_x_genes.genes.names.txt", cells = "/Users/kodeebao/Documents/StephaneMatrix/matrices/Apr2025_cold_rna_velocity/counts_unfiltered/cells_x_genes.barcodes.txt", feature.column = 1, mtx.transpose = TRUE)
COLD_NSC <- CreateSeuratObject(counts = expression_matrix, project = "COLD_NSC", min.cells = 3, min.features = 200)
COLD_NSC[["percent.mt"]] <- PercentageFeatureSet(COLD_NSC, pattern = "^mt-")
COLD_NSC <- subset(COLD_NSC, subset = nFeature_RNA > 1000 & percent.mt < 5)

expression_matrix <- ReadMtx(mtx = "/Users/kodeebao/Documents/StephaneMatrix/matrices/Apr2025_cold_rna_velocity/counts_unfiltered/cells_x_genes.ambiguous.mtx", features = "/Users/kodeebao/Documents/StephaneMatrix/matrices/Apr2025_cold_rna_velocity/counts_unfiltered/cells_x_genes.genes.names.txt", cells = "/Users/kodeebao/Documents/StephaneMatrix/matrices/Apr2025_cold_rna_velocity/counts_unfiltered/cells_x_genes.barcodes.txt", feature.column = 1, mtx.transpose = TRUE)
ambiguous_assay <- CreateAssayObject(counts = expression_matrix)
COLD_NSC[["ambiguous"]] <- subset(ambiguous_assay, cells=colnames(COLD_NSC))

expression_matrix <- ReadMtx(mtx = "/Users/kodeebao/Documents/StephaneMatrix/matrices/Apr2025_cold_rna_velocity/counts_unfiltered/cells_x_genes.mature.mtx", features = "/Users/kodeebao/Documents/StephaneMatrix/matrices/Apr2025_cold_rna_velocity/counts_unfiltered/cells_x_genes.genes.names.txt", cells = "/Users/kodeebao/Documents/StephaneMatrix/matrices/Apr2025_cold_rna_velocity/counts_unfiltered/cells_x_genes.barcodes.txt", feature.column = 1, mtx.transpose = TRUE)
spliced_assay <- CreateAssayObject(counts = expression_matrix)
COLD_NSC[["spliced"]] <- subset(spliced_assay, cells = colnames(COLD_NSC))

expression_matrix <- ReadMtx(mtx = "/Users/kodeebao/Documents/StephaneMatrix/matrices/Apr2025_cold_rna_velocity/counts_unfiltered/cells_x_genes.nascent.mtx", features = "/Users/kodeebao/Documents/StephaneMatrix/matrices/Apr2025_cold_rna_velocity/counts_unfiltered/cells_x_genes.genes.names.txt", cells = "/Users/kodeebao/Documents/StephaneMatrix/matrices/Apr2025_cold_rna_velocity/counts_unfiltered/cells_x_genes.barcodes.txt", feature.column = 1, mtx.transpose = TRUE)
unspliced_assay <- CreateAssayObject(counts = expression_matrix)
COLD_NSC[["unspliced"]] <- subset(unspliced_assay, cells = colnames(COLD_NSC))

# Set to 2 GiB
options(future.globals.maxSize = 2 * 1024^3)  

# SCT Transform (done instead of normalization, SF says it's better)
COLD_NSC <- SCTransform(object = COLD_NSC, assay = "RNA")
COLD_NSC <- SCTransform(object = COLD_NSC, assay = "ambiguous",  new.assay.name = "SCTambiguous")
COLD_NSC <- SCTransform(object = COLD_NSC, assay = "spliced",  new.assay.name = "SCTspliced")
COLD_NSC <- SCTransform(object = COLD_NSC, assay = "unspliced",  new.assay.name = "SCTunspliced")

# Set default assay
DefaultAssay(COLD_NSC) <- "SCT"
COLD_NSC <- RunPCA(COLD_NSC)

# number of cells
n_cells <- ncol(COLD_NSC)
print(n_cells)
# 12564 cells

saveRDS(COLD_NSC, "/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/COLD_NSC.rds")

# clear global environment before proceeding

# CON Cells #2 (Male and Female) ------------------------------------------

expression_matrix <- ReadMtx(mtx = "/Users/kodeebao/Documents/StephaneMatrix/matrices/May2025_control_2_rna_velocity/counts_unfiltered/cells_x_genes.cell.mtx", features = "/Users/kodeebao/Documents/StephaneMatrix/matrices/May2025_control_2_rna_velocity/counts_unfiltered/cells_x_genes.genes.names.txt", cells = "/Users/kodeebao/Documents/StephaneMatrix/matrices/May2025_control_2_rna_velocity/counts_unfiltered/cells_x_genes.barcodes.txt", feature.column = 1, mtx.transpose = TRUE)
CON_NSC <- CreateSeuratObject(counts = expression_matrix, project = "CON_NSC", min.cells = 3, min.features = 200)
CON_NSC[["percent.mt"]] <- PercentageFeatureSet(CON_NSC, pattern = "^mt-")
CON_NSC <- subset(CON_NSC, subset = nFeature_RNA > 1000 & percent.mt < 5)

expression_matrix <- ReadMtx(mtx = "/Users/kodeebao/Documents/StephaneMatrix/matrices/May2025_control_2_rna_velocity/counts_unfiltered/cells_x_genes.ambiguous.mtx", features = "/Users/kodeebao/Documents/StephaneMatrix/matrices/May2025_control_2_rna_velocity/counts_unfiltered/cells_x_genes.genes.names.txt", cells = "/Users/kodeebao/Documents/StephaneMatrix/matrices/May2025_control_2_rna_velocity/counts_unfiltered/cells_x_genes.barcodes.txt", feature.column = 1, mtx.transpose = TRUE)
ambiguous_assay <- CreateAssayObject(counts = expression_matrix)
CON_NSC[["ambiguous"]] <- subset(ambiguous_assay, cells=colnames(CON_NSC))

expression_matrix <- ReadMtx(mtx = "/Users/kodeebao/Documents/StephaneMatrix/matrices/May2025_control_2_rna_velocity/counts_unfiltered/cells_x_genes.mature.mtx", features = "/Users/kodeebao/Documents/StephaneMatrix/matrices/May2025_control_2_rna_velocity/counts_unfiltered/cells_x_genes.genes.names.txt", cells = "/Users/kodeebao/Documents/StephaneMatrix/matrices/May2025_control_2_rna_velocity/counts_unfiltered/cells_x_genes.barcodes.txt", feature.column = 1, mtx.transpose = TRUE)
spliced_assay <- CreateAssayObject(counts = expression_matrix)
CON_NSC[["spliced"]] <- subset(spliced_assay, cells = colnames(CON_NSC))

expression_matrix <- ReadMtx(mtx = "/Users/kodeebao/Documents/StephaneMatrix/matrices/May2025_control_2_rna_velocity/counts_unfiltered/cells_x_genes.nascent.mtx", features = "/Users/kodeebao/Documents/StephaneMatrix/matrices/May2025_control_2_rna_velocity/counts_unfiltered/cells_x_genes.genes.names.txt", cells = "/Users/kodeebao/Documents/StephaneMatrix/matrices/May2025_control_2_rna_velocity/counts_unfiltered/cells_x_genes.barcodes.txt", feature.column = 1, mtx.transpose = TRUE)
unspliced_assay <- CreateAssayObject(counts = expression_matrix)
CON_NSC[["unspliced"]] <- subset(unspliced_assay, cells = colnames(CON_NSC))

# Set to 2 GiB
options(future.globals.maxSize = 2 * 1024^3)  

# SCT Transform (done instead of normalization, SF says it's better)
CON_NSC <- SCTransform(object = CON_NSC, assay = "RNA")
CON_NSC <- SCTransform(object = CON_NSC, assay = "ambiguous",  new.assay.name = "SCTambiguous")
CON_NSC <- SCTransform(object = CON_NSC, assay = "spliced",  new.assay.name = "SCTspliced")
CON_NSC <- SCTransform(object = CON_NSC, assay = "unspliced",  new.assay.name = "SCTunspliced")

# Set default assay
DefaultAssay(CON_NSC) <- "SCT"
CON_NSC <- RunPCA(CON_NSC)

# number of cells
n_cells <- ncol(CON_NSC)
print(n_cells)
# 10079 cells

saveRDS(CON_NSC, "/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/CON_NSC.rds")

# clear global environment before proceeding

# merging objects  -------------------------------------------------
COLD_NSC <- readRDS('/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/COLD_NSC.rds')
CON_NSC <- readRDS('/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/CON_NSC.rds')

# could not merge so trying DietSeurat, exclude the SCTambiguous, SCTspliced, SCTunspliced assays from both seurat objects I want to merge, then I merge them
# will need to check if SF needs them post-clustered before doing that so he can run RNA Velocity, SF said it is fine and will still be able to get results 

# Define the assays to exclude
exclude_assays <- c("SCTambiguous", "SCTspliced", "SCTunspliced")

# Get the assays to keep for each object
keep_COLD_assays <- setdiff(Assays(COLD_NSC), exclude_assays)
keep_CON_assays <- setdiff(Assays(CON_NSC), exclude_assays)

# Slim down the Seurat objects
COLD_NSC_diet <- DietSeurat(COLD_NSC, assays = keep_COLD_assays, misc = TRUE)
CON_NSC_diet <- DietSeurat(CON_NSC, assays = keep_CON_assays, misc = TRUE)

# save diet objects
saveRDS(COLD_NSC_diet, "/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/COLD_NSC_diet.rds")
saveRDS(CON_NSC_diet, "/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/CON_NSC_diet.rds")

# clear global environment before proceeding 

# load in diet objects 
COLD_NSC_diet <- readRDS('/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/COLD_NSC_diet.rds')
CON_NSC_diet <- readRDS('/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/CON_NSC_diet.rds')

# Merge COLD_NSC and CON_NSC
merged_object <- merge(
  x = CON_NSC_diet,                          # reference object
  y = COLD_NSC_diet,                         # object to merge
  add.cell.ids = c("CON", "COLD"),     # prefixes to distinguish cells
  project = "Merged_NSC")

saveRDS(merged_object, "/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/merged_object.rds")

# clear global environment before proceeding 

# analysis without integration --------------------------------------------

# load in merged_object
merged_object <- readRDS('/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/merged_object.rds')

# Set to 2 GiB
options(future.globals.maxSize = 2 * 1024^3)  

# SCTransform merged object
merged_object <- SCTransform(object = merged_object, assay = "RNA")
merged_object <- SCTransform(object = merged_object, assay = "ambiguous",  new.assay.name = "SCTambiguous")
merged_object <- SCTransform(object = merged_object, assay = "spliced",  new.assay.name = "SCTspliced")
merged_object <- SCTransform(object = merged_object, assay = "unspliced",  new.assay.name = "SCTunspliced")

DefaultAssay(merged_object) <- "SCT"
merged_object <- RunPCA(merged_object)
ElbowPlot(merged_object, ndims = 40, reduction = "pca")		

saveRDS(merged_object, "/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/merged_object_SCTransformed.rds")


# clustering merged object that is SCTransformed --------------------------
# merged object contains data sets COLD M & F and CON M & F, they have been dieted then SCTransformed again

# load in merged_object
merged_object <- readRDS('/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/merged_object_SCTransformed.rds')

# Confirm number of cells
ncol(merged_object)
# 22643 cells

# now JR wants to pull out Vim+ cells and make clustrees
# cluster then pull out Vim+ cells because you want to cluster them based on all the cells first 

# pull out Vim+ cells from merged object ----------------------------------

# Subset cells expressing vimentin (Vim > 0)
merged_object_vim_pos <- subset(merged_object, subset = Vim > 0)

# Confirm number of cells
ncol(merged_object_vim_pos)
# 13624 cells

DefaultAssay(merged_object_vim_pos) <- "SCT"
merged_object_vim_pos <- RunPCA(merged_object_vim_pos)
ElbowPlot(merged_object_vim_pos, ndims = 40, reduction = "pca")		

saveRDS(merged_object_vim_pos, "/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/merged_object_vim_pos.rds")


# load in Vim+ subsetted cells --------------------------------------------

merged_object_vim_pos <- readRDS('/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/merged_object_vim_pos.rds')

# Confirm number of cells
ncol(merged_object_vim_pos)
# 13624 cells

# using 35 dims, 0.1 res to cluster (found on clustree on separate workflow "Testing Parameters to cluster COLD and CON datasets v2, Vim+ only.R" and discussed with JR)
# 0.1 res looked to have more cells in cluster 1 (fewer red cells in cluster 1)
merged_object_vim_pos <- FindNeighbors(merged_object_vim_pos, dims = 1:35)	
merged_object_vim_pos <- FindClusters(merged_object_vim_pos, resolution = 0.1)		
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
#   COLD_NSC     CON_NSC 
#    3128        2543      (Total 5671)

DefaultAssay(merged_object_vim_pos_subset) <- "SCT"
merged_object_vim_pos_subset <- RunPCA(merged_object_vim_pos_subset)
ElbowPlot(merged_object_vim_pos_subset, ndims = 40, reduction = "pca")		

saveRDS(merged_object_vim_pos_subset, "/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/merged_object_vim_pos_subset.rds")

# load in Vim+ cells that have been subsetted (cluster 1 or Ascl1+) -----------------------------------------

merged_object_vim_pos_subset <- readRDS('/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/merged_object_vim_pos_subset.rds')

# testing parameters on cluster 1 + Ascl1+ subset on separate workflow 

# using 30 dims, 0.2 res -- confirmed with JR -----------------------------
merged_object_vim_pos_subset <- FindNeighbors(merged_object_vim_pos_subset, dims = 1:30)	
merged_object_vim_pos_subset <- FindClusters(merged_object_vim_pos_subset, resolution = 0.2)		
merged_object_vim_pos_subset <- RunUMAP(merged_object_vim_pos_subset, dims = 1:30)
DimPlot(merged_object_vim_pos_subset, reduction = "umap")
DimPlot(merged_object_vim_pos_subset, reduction = "umap", group.by = c("orig.ident", "seurat_clusters"), split.by = "orig.ident")

# quality check  ----------------------------------------------------------

# Total number of cells 
ncol(merged_object_vim_pos_subset)
# 5761

# Number of cells in each cluster 
table(Idents(merged_object_vim_pos_subset))
#    0    1    2    3    4    5 
# 2275 1478  954  525  342   97 

# Number of CON vs. COLD cells in total
table(merged_object_vim_pos_subset$orig.ident)
#    COLD_NSC    CON_NSC 
#    3128        2543 

# Number of CON and COLD cells in each cluster
table(merged_object_vim_pos_subset$seurat_clusters, merged_object_vim_pos_subset$orig.ident)
#      COLD_NSC     CON_NSC
# 0      1214       1061
# 1      760        718
# 2      665        289
# 3      232        293
# 4      218        124
# 5      39         58


# removing male and female double positive and non-expressing cells -------
# From "Analysis on Merged CON & COLD" workflow 
# using "Xist", "Kdm6a", "Eif2s3x" for females and "Ddx3y", "Uty" for males gave cells overlapping in both males and females
# females = COLD_NSC = 2103; CON_NSC = 2509 // males = COLD_NSC = 971; CON_NSC = 0
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
# 2364   2121    (Total 4485)

# Add descriptive cluster names (if using custom cluster identities)
clean_combined_object$cluster_names <- recode(
  as.character(clean_combined_object$seurat_clusters),
  "0" = "Neuronal Precursors (Ascl1+)",
  "1" = "Glial Precursors",
  "2" = "Cycling NSPCs",
  "3" = "NPCs-Glia",
  "4" = "NPCs-Neuronal",
  "5" = "Neuronal Precursors (Neurog2+)")

# Define custom cluster colors
cluster_colors <- c(
  "Neuronal Precursors (Ascl1+)" = "#F8766D",
  "Glial Precursors" = "#B79F00",
  "Cycling NSPCs" = "#00BA38",
  "NPCs-Glia" = "#00BFC4",
  "NPCs-Neuronal" = "#649DFF",
  "Neuronal Precursors (Neurog2+)" = "#F564E3")

# store treatment in metadata
clean_combined_object$treatment <- ifelse(grepl("COLD", clean_combined_object$orig.ident), "COLD", "CON")

# Plot UMAPs using cluster_names and custom colors
DimPlot(clean_combined_object, reduction = "umap", group.by = "cluster_names", cols = cluster_colors) + ggtitle("Combined M and F - CON and COLD")

table(clean_combined_object$sex)
# Female   Male 
# 2364    2121 

table(clean_combined_object$treatment)
# COLD  CON 
# 2392 2093 

table(clean_combined_object$cluster_names)
# Cycling NSPCs               Glial Precursors   Neuronal Precursors (Ascl1+) 
#  756                           1137                           1839 
# Neuronal Precursors (Neurog2+)                      NPCs-Glia                  NPCs-Neuronal 
#              77                                      408                            268 

saveRDS(clean_combined_object, file = "/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/clean_combined_object_CON&COLD_strat2.rds")


# load in clean combined object -------------------------------------------

clean_combined_object <- readRDS('/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/clean_combined_object_CON&COLD_strat2.rds')


# splitting "clean_combined_object" into male and female objects ----------

# Subset cells
female_obj <- subset(clean_combined_object, subset = sex == "Female")
male_obj <- subset(clean_combined_object, subset = sex == "Male")

# check
cat("Female cells:", ncol(female_obj), "\n") # 2364
cat("Male cells:", ncol(male_obj), "\n") # 2121

# save 
saveRDS(female_obj, file = "/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/clean_female_obj_CON&COLD_strat2.rds")
saveRDS(male_obj, file = "/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/clean_male_obj_CON&COLD_strat2.rds")

# load 
female_obj <- readRDS('/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/clean_female_obj_CON&COLD_strat2.rds')
male_obj <- readRDS('/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/clean_male_obj_CON&COLD_strat2.rds')


# splitting "clean_combined_object" into CON and COLD objects (has both male and female) ----------

# Subset cells
con_obj <- subset(clean_combined_object, subset = treatment == "CON")
cold_obj <- subset(clean_combined_object, subset = treatment == "COLD")

# check
cat("CON cells:", ncol(con_obj), "\n") # 2093
cat("COLD cells:", ncol(cold_obj), "\n") # 2392

# save 
saveRDS(con_obj, file = "/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/clean_con_obj_M&F_strat2.rds")
saveRDS(cold_obj, file = "/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/clean_cold_obj_M&F_strat2.rds")

# load 
con_obj <- readRDS('/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/clean_con_obj_M&F_strat2.rds')
cold_obj <- readRDS('/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/clean_cold_obj_M&F_strat2.rds')


# splitting "clean_combined_object" into by sex and CON/COLD ----------

# Female + CON
female_con_obj <- subset(clean_combined_object, subset = sex == "Female" & treatment == "CON")

# Female + COLD
female_cold_obj <- subset(clean_combined_object, subset = sex == "Female" & treatment == "COLD")

# Male + CON
male_con_obj <- subset(clean_combined_object, subset = sex == "Male" & treatment == "CON")

# Male + COLD
male_cold_obj <- subset(clean_combined_object, subset = sex == "Male" & treatment == "COLD")

# check
sapply(list(female_con_obj, female_cold_obj, male_con_obj, male_cold_obj), ncol) # 1167, 1197, 926, 1195

# save
saveRDS(female_con_obj, file = "/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/female_con_obj.rds")
saveRDS(female_cold_obj, file = "/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/female_cold_obj.rds")
saveRDS(male_con_obj, file = "/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/male_con_obj.rds")
saveRDS(male_cold_obj, file = "/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/male_cold_obj.rds")

# load
female_con_obj <- readRDS('/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/female_con_obj.rds')
female_cold_obj <- readRDS('/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/female_cold_obj.rds')
male_con_obj <- readRDS('/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/male_con_obj.rds')
male_cold_obj <- readRDS('/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/male_cold_obj.rds')


# make umaps --------------------------------------------------------------

# Generate UMAPs using DimPlot
umap_female_con <- DimPlot(female_con_obj, reduction = "umap", group.by = "seurat_clusters") + 
  ggtitle("Female - Control")

umap_female_cold <- DimPlot(female_cold_obj, reduction = "umap", group.by = "seurat_clusters") + 
  ggtitle("Female - Cold")

umap_male_con <- DimPlot(male_con_obj, reduction = "umap", group.by = "seurat_clusters") + 
  ggtitle("Male - Control")

umap_male_cold <- DimPlot(male_cold_obj, reduction = "umap", group.by = "seurat_clusters") + 
  ggtitle("Male - Cold")

# Combine UMAPs into a single panel
umap_female_con | umap_female_cold | umap_male_con | umap_male_cold
