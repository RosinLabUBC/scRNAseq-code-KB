# going to try and remove all the CON F final cells I had from my original dataset from the raw data then merge with the CON F dataset verification and process as usual 


library(Seurat)
library(dplyr)
library(patchwork)
library(clustree)


# original CON cells ------------------------------------------------------
# load in CON dataset from original processing
CON_NSC <- readRDS('/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/CON_NSC.rds')

# number of cells
n_cells <- ncol(CON_NSC)
print(n_cells)
# 10079 cells

# remove CON F cells
female_con_obj <- readRDS('/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/female_con_obj.rds')

# number of cells
n_cells <- ncol(female_con_obj)
print(n_cells)
# 1167 cells

# Remove "CON_" prefix from female_con_obj cell names
female_con_obj <- RenameCells(female_con_obj, new.names = sub("^CON_", "", Cells(female_con_obj)))

# Get the cell names to remove
cells_to_remove <- Cells(female_con_obj)

# Subset CON_NSC to exclude those cells
CON_NSC_no_CON_F <- subset(
  CON_NSC,
  cells = setdiff(Cells(CON_NSC), cells_to_remove))

# Check number of cells
length(Cells(CON_NSC))          # original = 10079
length(Cells(CON_NSC_no_CON_F)) # after removal = 8912

saveRDS(CON_NSC_no_CON_F, "/Users/kodeebao/Documents/Strategy 2 Final Analysis/CON F Verification 2/saved objects/CON_NSC_original_no_CON_F.rds")

# Process CON Cells #1 (Female ONLY) --------------------------------------------

expression_matrix <- ReadMtx(mtx = "/Users/kodeebao/Documents/StephaneMatrix/matrices/Mar2025_control_rna_velocity/counts_unfiltered/cells_x_genes.cell.mtx", features = "/Users/kodeebao/Documents/StephaneMatrix/matrices/Mar2025_control_rna_velocity/counts_unfiltered/cells_x_genes.genes.names.txt", cells = "/Users/kodeebao/Documents/StephaneMatrix/matrices/Mar2025_control_rna_velocity/counts_unfiltered/cells_x_genes.barcodes.txt", feature.column = 1, mtx.transpose = TRUE)
CON_NSC <- CreateSeuratObject(counts = expression_matrix, project = "CON_NSC", min.cells = 3, min.features = 200)
CON_NSC[["percent.mt"]] <- PercentageFeatureSet(CON_NSC, pattern = "^mt-")
CON_NSC <- subset(CON_NSC, subset = nFeature_RNA > 1000 & percent.mt < 5)

expression_matrix <- ReadMtx(mtx = "/Users/kodeebao/Documents/StephaneMatrix/matrices/Mar2025_control_rna_velocity/counts_unfiltered/cells_x_genes.ambiguous.mtx", features = "/Users/kodeebao/Documents/StephaneMatrix/matrices/Mar2025_control_rna_velocity/counts_unfiltered/cells_x_genes.genes.names.txt", cells = "/Users/kodeebao/Documents/StephaneMatrix/matrices/Mar2025_control_rna_velocity/counts_unfiltered/cells_x_genes.barcodes.txt", feature.column = 1, mtx.transpose = TRUE)
ambiguous_assay <- CreateAssayObject(counts = expression_matrix)
CON_NSC[["ambiguous"]] <- subset(ambiguous_assay, cells=colnames(CON_NSC))

expression_matrix <- ReadMtx(mtx = "/Users/kodeebao/Documents/StephaneMatrix/matrices/Mar2025_control_rna_velocity/counts_unfiltered/cells_x_genes.mature.mtx", features = "/Users/kodeebao/Documents/StephaneMatrix/matrices/Mar2025_control_rna_velocity/counts_unfiltered/cells_x_genes.genes.names.txt", cells = "/Users/kodeebao/Documents/StephaneMatrix/matrices/Mar2025_control_rna_velocity/counts_unfiltered/cells_x_genes.barcodes.txt", feature.column = 1, mtx.transpose = TRUE)
spliced_assay <- CreateAssayObject(counts = expression_matrix)
CON_NSC[["spliced"]] <- subset(spliced_assay, cells = colnames(CON_NSC))

expression_matrix <- ReadMtx(mtx = "/Users/kodeebao/Documents/StephaneMatrix/matrices/Mar2025_control_rna_velocity/counts_unfiltered/cells_x_genes.nascent.mtx", features = "/Users/kodeebao/Documents/StephaneMatrix/matrices/Mar2025_control_rna_velocity/counts_unfiltered/cells_x_genes.genes.names.txt", cells = "/Users/kodeebao/Documents/StephaneMatrix/matrices/Mar2025_control_rna_velocity/counts_unfiltered/cells_x_genes.barcodes.txt", feature.column = 1, mtx.transpose = TRUE)
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
# 13455 cells

saveRDS(CON_NSC, "/Users/kodeebao/Documents/Strategy 2 Final Analysis/CON F Verification 2/saved objects/CON_NSC_female_verification_F_only.rds")

# clear global environment before proceeding


# load objects, diet, then merge ------------------------------------------------------------
CON_NSC_orig_no_F <- readRDS('/Users/kodeebao/Documents/Strategy 2 Final Analysis/CON F Verification 2/saved objects/CON_NSC_original_no_CON_F.rds')

CON_NSC_F_only <- readRDS('/Users/kodeebao/Documents/Strategy 2 Final Analysis/CON F Verification 2/saved objects/CON_NSC_female_verification_F_only.rds')

# Define the assays to exclude
exclude_assays <- c("SCTambiguous", "SCTspliced", "SCTunspliced")

# Get the assays to keep for CON object
keep_CON_assays_orig <- setdiff(Assays(CON_NSC_orig_no_F), exclude_assays)
keep_CON_assays_verif <- setdiff(Assays(CON_NSC_F_only), exclude_assays)

# Slim down the Seurat objects
CON_NSC_orig_no_F_diet <- DietSeurat(CON_NSC_orig_no_F, assays = keep_CON_assays_orig, misc = TRUE)
CON_NSC_F_only_diet <- DietSeurat(CON_NSC_F_only, assays = keep_CON_assays_verif, misc = TRUE)

# save diet objects
saveRDS(CON_NSC_orig_no_F_diet, "/Users/kodeebao/Documents/Strategy 2 Final Analysis/CON F Verification 2/saved objects/CON_NSC_orig_diet.rds")
saveRDS(CON_NSC_F_only_diet, "/Users/kodeebao/Documents/Strategy 2 Final Analysis/CON F Verification 2/saved objects/CON_NSC_verif_diet.rds")


# clear global environment before proceeding 

# load in diet objects 
COLD_NSC_diet <- readRDS('/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/COLD_NSC_diet.rds')
CON_NSC_orig_diet <- readRDS('/Users/kodeebao/Documents/Strategy 2 Final Analysis/CON F Verification 2/saved objects/CON_NSC_orig_diet.rds')
CON_NSC_verif_diet <- readRDS('/Users/kodeebao/Documents/Strategy 2 Final Analysis/CON F Verification 2/saved objects/CON_NSC_verif_diet.rds')

# Merge the three Seurat objects with unique prefixes
merged_object <- merge(
  x = CON_NSC_orig_diet,                              # reference object
  y = list(CON_NSC_verif_diet, COLD_NSC_diet),        # objects to merge
  add.cell.ids = c("CON_orig", "CON_verif", "COLD"),  # unique prefixes
  project = "Merged_NSC")

# Create a 'treatment' metadata column: group all CON cells together
merged_object$treatment <- ifelse(
  grepl("^CON", Cells(merged_object)),  # any cell starting with "CON"
  "CON", "COLD")

saveRDS(merged_object, "/Users/kodeebao/Documents/Strategy 2 Final Analysis/CON F Verification 2/saved objects/merged_object.rds")

# Using ncol on the Seurat object
n_cells <- ncol(merged_object)
print(n_cells)
# 34931 cells

# clear global environment before proceeding 

# analysis without integration --------------------------------------------

# load in merged_object
merged_object <- readRDS('/Users/kodeebao/Documents/Strategy 2 Final Analysis/CON F Verification 2/saved objects/merged_object.rds')

# Set to 2 GiB
options(future.globals.maxSize = 2 * 1024^3)  

# SCTransform merged object
merged_object <- SCTransform(object = merged_object, assay = "RNA")
# merged_object <- SCTransform(object = merged_object, assay = "ambiguous",  new.assay.name = "SCTambiguous")
# merged_object <- SCTransform(object = merged_object, assay = "spliced",  new.assay.name = "SCTspliced")
# merged_object <- SCTransform(object = merged_object, assay = "unspliced",  new.assay.name = "SCTunspliced")
# don't have memory to SCTransform ambiguous, spliced, or unspliced, but don't need because will just use RNA -- just to verify DEG which uses RNA assay, the others are for RNA velocity

DefaultAssay(merged_object) <- "SCT"
merged_object <- RunPCA(merged_object)
ElbowPlot(merged_object, ndims = 40, reduction = "pca")		

saveRDS(merged_object, "/Users/kodeebao/Documents/Strategy 2 Final Analysis/CON F Verification 2/saved objects/merged_object_SCTransformed.rds")

# clustering merged object that is SCTransformed --------------------------
# merged object contains data sets COLD M & F and CON F (from dataset 1) and CON cells from dataset 1 (without any of the final female cells from original), they have been dieted then SCTransformed again

# load in merged_object
merged_object <- readRDS('/Users/kodeebao/Documents/Strategy 2 Final Analysis/CON F Verification 2/saved objects/merged_object_SCTransformed.rds')

# Confirm number of cells
ncol(merged_object)
# 34931 cells

# pull out Vim+ cells from merged object ----------------------------------

# Subset cells expressing vimentin (Vim > 0)
merged_object_vim_pos <- subset(merged_object, subset = Vim > 0)

# Confirm number of cells
ncol(merged_object_vim_pos)
# 21132 cells

DefaultAssay(merged_object_vim_pos) <- "SCT"
merged_object_vim_pos <- RunPCA(merged_object_vim_pos)
ElbowPlot(merged_object_vim_pos, ndims = 40, reduction = "pca")		

saveRDS(merged_object_vim_pos, "/Users/kodeebao/Documents/Strategy 2 Final Analysis/CON F Verification 2/saved objects/merged_object_vim_pos.rds")

# load in Vim+ subsetted cells --------------------------------------------

merged_object_vim_pos <- readRDS('/Users/kodeebao/Documents/Strategy 2 Final Analysis/CON F Verification 2/saved objects/merged_object_vim_pos.rds')

# Confirm number of cells
ncol(merged_object_vim_pos)
# 21132 cells

# using 35 dims, 0.1 res to cluster (was used on regular workflow 'workflow to cluster and subset.R')
merged_object_vim_pos <- FindNeighbors(merged_object_vim_pos, dims = 1:35)	
merged_object_vim_pos <- FindClusters(merged_object_vim_pos, resolution = 0.1)		
merged_object_vim_pos <- RunUMAP(merged_object_vim_pos, dims = 1:35)
DimPlot(merged_object_vim_pos, reduction = "umap")
DimPlot(merged_object_vim_pos, reduction = "umap", group.by = c("orig.ident", "seurat_clusters"), split.by = "orig.ident")

# looks very similar to the original workflow, so based on testing gene markers, going to subset out cluster 1 and all Ascl1+ cells

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
#    3125        4697      (Total 7822)

DefaultAssay(merged_object_vim_pos_subset) <- "SCT"
merged_object_vim_pos_subset <- RunPCA(merged_object_vim_pos_subset)
ElbowPlot(merged_object_vim_pos_subset, ndims = 40, reduction = "pca")		

saveRDS(merged_object_vim_pos_subset, "/Users/kodeebao/Documents/Strategy 2 Final Analysis/CON F Verification 2/saved objects/merged_object_vim_pos_subset.rds")

# load in Vim+ cells that have been subsetted (cluster 1 or Ascl1+) -----------------------------------------

merged_object_vim_pos_subset <- readRDS('/Users/kodeebao/Documents/Strategy 2 Final Analysis/CON F Verification 2/saved objects/merged_object_vim_pos_subset.rds')

# using 30 dims, 0.15 res (30 dims, 0.2 res was used on regular workflow 'workflow to cluster and subset.R') ----------------
merged_object_vim_pos_subset <- FindNeighbors(merged_object_vim_pos_subset, dims = 1:30)	
merged_object_vim_pos_subset <- FindClusters(merged_object_vim_pos_subset, resolution = 0.18)		
merged_object_vim_pos_subset <- RunUMAP(merged_object_vim_pos_subset, dims = 1:30)
DimPlot(merged_object_vim_pos_subset, reduction = "umap")
DimPlot(merged_object_vim_pos_subset, reduction = "umap", group.by = c("orig.ident", "seurat_clusters"), split.by = "orig.ident")

# using 30 dims, 18 res because larger set of cells made clustering a bit different for cycling NSPCs (instead of 30 dims, 0.2 res) and looks exactly how we want it to

# cluster 2 is Cycling NSPCs, based on markers (can be found under 'testing gene markers 2.R' workflow under CON F Verification 2)

# removing male and female double positive and non-expressing cells -------
# like in regular workflow

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
# 4470    2145    (Total 6615)

# add cluster 6 into cluster 3, based on markers represent NPCs-Glia population
clean_combined_object <- RenameIdents(
  clean_combined_object,
  `6` = "3")

clean_combined_object$cluster_names <- as.character(
  Idents(clean_combined_object))

# Add descriptive cluster names (if using custom cluster identities)
clean_combined_object$cluster_names <- dplyr::recode(
  clean_combined_object$cluster_names,
  "2" = "Cycling NSPCs",
  .default = clean_combined_object$cluster_names)

# Define custom cluster colors
cluster_colors1 <- c(
  "Cycling NSPCs" = "#00BA38")

# Define custom cluster colors
cluster_colors <- c(
  "2" = "#00BA38",
  "0" = "#F8766D",
  "1" = "#B79F00",
  "3" = "#00BFC4",
  "4" = "#649DFF",
  "5" = "#F564E3")


# store treatment in metadata
# clean_combined_object$treatment <- ifelse(grepl("COLD", clean_combined_object$orig.ident), "COLD", "CON") 
# already stored previously

DimPlot(clean_combined_object, reduction = "umap", cols = cluster_colors)

# Plot UMAPs using cluster_names and custom colors
DimPlot(clean_combined_object, reduction = "umap", group.by = "cluster_names", cols = cluster_colors1) + ggtitle("Combined M and F - CON and COLD")

table(clean_combined_object$sex)
# Female   Male 
# 4470    2145 

table(clean_combined_object$treatment)
# COLD  CON 
# 2391 4224 

table(clean_combined_object$cluster_names)
# Cycling NSPCs - 1188

saveRDS(clean_combined_object, file = "/Users/kodeebao/Documents/Strategy 2 Final Analysis/CON F Verification 2/saved objects/clean_combined_object_CON&COLD_F_verification.rds")

# load in clean combined object -------------------------------------------

clean_combined_object <- readRDS('/Users/kodeebao/Documents/Strategy 2 Final Analysis/CON F Verification 2/saved objects/clean_combined_object_CON&COLD_F_verification.rds')

# check "clean_combined_object" into by sex and CON/COLD ----------

# Female + CON
female_con_obj <- subset(clean_combined_object, subset = sex == "Female" & treatment == "CON")

# Female + COLD
female_cold_obj <- subset(clean_combined_object, subset = sex == "Female" & treatment == "COLD")

# Male + CON
male_con_obj <- subset(clean_combined_object, subset = sex == "Male" & treatment == "CON")

# Male + COLD
male_cold_obj <- subset(clean_combined_object, subset = sex == "Male" & treatment == "COLD")

# check
sapply(list(female_con_obj, female_cold_obj, male_con_obj, male_cold_obj), ncol) # 3274, 1196, 950, 1195

# save
saveRDS(female_con_obj, file = "/Users/kodeebao/Documents/Strategy 2 Final Analysis/CON F Verification 2/saved objects/female_con_obj.rds")
saveRDS(female_cold_obj, file = "/Users/kodeebao/Documents/Strategy 2 Final Analysis/CON F Verification 2/saved objects/female_cold_obj.rds")
saveRDS(male_con_obj, file = "/Users/kodeebao/Documents/Strategy 2 Final Analysis/CON F Verification 2/saved objects/male_con_obj.rds")
saveRDS(male_cold_obj, file = "/Users/kodeebao/Documents/Strategy 2 Final Analysis/CON F Verification 2/saved objects/male_cold_obj.rds")


# make umaps --------------------------------------------------------------

# Generate UMAPs using DimPlot
umap_female_con <- DimPlot(female_con_obj, reduction = "umap", cols = cluster_colors) + 
  ggtitle("Female - Control")

umap_female_cold <- DimPlot(female_cold_obj, reduction = "umap", cols = cluster_colors) + 
  ggtitle("Female - Cold")

umap_male_con <- DimPlot(male_con_obj, reduction = "umap", cols = cluster_colors) + 
  ggtitle("Male - Control")

umap_male_cold <- DimPlot(male_cold_obj, reduction = "umap", cols = cluster_colors) + 
  ggtitle("Male - Cold")

# Combine UMAPs into a single panel
umap_female_con | umap_female_cold | umap_male_con | umap_male_cold

# Birc5
plot <- FeaturePlot(clean_combined_object, features = c("Birc5"))
ggsave(filename = '/Users/kodeebao/Documents/Strategy 2 Final Analysis/CON F Verification 2/Analysis on Vim+ subsetted cells/umaps/Birc5_vim_pos_subset_feature.png', plot, width = 8, height = 5)

plot <- VlnPlot(clean_combined_object, features = c("Birc5"))
ggsave(filename = '/Users/kodeebao/Documents/Strategy 2 Final Analysis/CON F Verification/Analysis on Vim+ subsetted cells/umaps/Birc5_vim_pos_subset_vln.png', plot, width = 8, height = 5)

# cycling markers
plot <- FeaturePlot(clean_combined_object, features = c("Birc5", 
                                                        "Ccnb1",
                                                        "Sgo2a",
                                                        "Cdca8",
                                                        "Tpx2",
                                                        "Ube2c"))
ggsave(filename = '/Users/kodeebao/Documents/Strategy 2 Final Analysis/CON F Verification 2/Analysis on Vim+ subsetted cells/umaps/cycling_markers_feature.png', plot, width = 10, height = 12)


# figure out cluster 6  -----------------------------------------------------------

# npcs-glia and npcs-neuronal cell markers
plot <- FeaturePlot(clean_combined_object, features = c("Olig1",
                                                        "Olig2",
                                                        "Pdgfra",
                                                        "Dlx1",
                                                        "Dlx2",
                                                        "Arx"))
ggsave(filename = '/Users/kodeebao/Documents/Strategy 2 Final Analysis/CON F Verification 2/Analysis on Vim+ subsetted cells/umaps/npcs-glia and neuronal_markers_feature.png', plot, width = 10, height = 10)

