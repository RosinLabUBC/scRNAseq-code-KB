
# analysis 1 --------------------------------------------------------------

# merging blackshaw with our CON and COLD data

library(Seurat)
library(dplyr)
library(patchwork)
library(clustree)

# Set to 2 GiB
options(future.globals.maxSize = 2 * 1024^3)  

# load in blackshaw objects 
Rep1_2025 <- readRDS('/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/saved objects/Rep1_2025.rds')
Rep2_2025 <- readRDS('/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/saved objects/Rep2_2025.rds')
Rep1_2020 <- readRDS('/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/saved objects/Rep1_2020.rds')

# load in KB diet objects 
COLD_NSC_diet <- readRDS('/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/COLD_NSC_diet.rds')
CON_NSC_diet <- readRDS('/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/CON_NSC_diet.rds')

# merge all objects
merged_object <- merge(
  x = CON_NSC_diet,                          # reference object
  y = list(COLD_NSC_diet, Rep1_2025, Rep2_2025, Rep1_2020),                         # object to merge
  add.cell.ids = c("CON_KB", "COLD_KB", "CON_blackshaw_1", "CON_blackshaw_2", "CON_blackshaw_3"),     # prefixes to distinguish cells
  project = "Merged_NSC")

saveRDS(merged_object, "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/blackshaw + kb/merged_object.rds")

# checked umap on merged object before SCTransforming and the cells were distributed nicely among eachother, don't need to integrate 

# analysis without integration --------------------------------------------
#SCTransform merged_object

# load in merged_object
merged_object <- readRDS('/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/blackshaw + kb/merged_object.rds')

# Set to 2 GiB
options(future.globals.maxSize = 2 * 1024^3)  

# SCTransform merged object
merged_object <- SCTransform(object = merged_object, assay = "RNA")

DefaultAssay(merged_object) <- "SCT"
merged_object <- RunPCA(merged_object)
ElbowPlot(merged_object, ndims = 40, reduction = "pca")		

saveRDS(merged_object, "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/blackshaw + kb/merged_object_SCTransformed.rds")

# clustering merged object that is SCTransformed --------------------------
# merged object contains KB CON and COLD data + 3 blackshaw CON data sets, they have been SCTransformed again

# load in merged_object
merged_object <- readRDS('/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/blackshaw + kb/merged_object_SCTransformed.rds')

# Confirm number of cells
ncol(merged_object)
# 34328 cells

# pull out Vim+ cells from merged object ----------------------------------

# Subset cells expressing vimentin (Vim > 0)
merged_object_vim_pos <- subset(merged_object, subset = Vim > 0)

# Confirm number of cells
ncol(merged_object_vim_pos)
# 19195 cells

DefaultAssay(merged_object_vim_pos) <- "SCT"
merged_object_vim_pos <- RunPCA(merged_object_vim_pos)
ElbowPlot(merged_object_vim_pos, ndims = 40, reduction = "pca")		

saveRDS(merged_object_vim_pos, "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/blackshaw + kb/merged_object_vim_pos.rds")

# load in Vim+ subsetted cells --------------------------------------------
merged_object_vim_pos <- readRDS('/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/blackshaw + kb/merged_object_vim_pos.rds')

# Confirm number of cells
ncol(merged_object_vim_pos)
# 19195 cells

# using 35 dims, 0.05 res to cluster 
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
#      2020_E15_Rep1   2025_E15.5_Rep1   2025_E15.5_Rep2        COLD_NSC         CON_NSC 
#            58            1136               1330               3127            2550 
# 8201 total cells

DefaultAssay(merged_object_vim_pos_subset) <- "SCT"
merged_object_vim_pos_subset <- RunPCA(merged_object_vim_pos_subset)
ElbowPlot(merged_object_vim_pos_subset, ndims = 40, reduction = "pca")		

saveRDS(merged_object_vim_pos_subset, "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/blackshaw + kb/merged_object_vim_pos_subset.rds")

# load in Vim+ cells that have been subsetted (cluster 1 or Ascl1+) -----------------------------------------

merged_object_vim_pos_subset <- readRDS('/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/blackshaw + kb/merged_object_vim_pos_subset.rds')

# Total number of cells 
ncol(merged_object_vim_pos_subset)
# 8201 cells

# tested parameters on cluster 1 + Ascl1+ subset on clustree

# using 30 dims, 0.1 res -----------------------------
merged_object_vim_pos_subset <- FindNeighbors(merged_object_vim_pos_subset, dims = 1:30)	
merged_object_vim_pos_subset <- FindClusters(merged_object_vim_pos_subset, resolution = 0.1)		
merged_object_vim_pos_subset <- RunUMAP(merged_object_vim_pos_subset, dims = 1:30)
DimPlot(merged_object_vim_pos_subset, reduction = "umap")
DimPlot(merged_object_vim_pos_subset, reduction = "umap", group.by = c("orig.ident", "seurat_clusters"), split.by = "orig.ident")

# quality check  ----------------------------------------------------------

# Total number of cells 
ncol(merged_object_vim_pos_subset)
# 8201

# Number of cells in each cluster 
table(Idents(merged_object_vim_pos_subset))
#    0     1     2     3     4    5 
#  2847  2658  1863  562   216   55 

# Number of CON vs. COLD cells in total
table(merged_object_vim_pos_subset$orig.ident)
#    2020_E15_Rep1    2025_E15.5_Rep1   2025_E15.5_Rep2       COLD_NSC      CON_NSC 
#        58            1136              1330                 3127            2550 

# Number of CON and COLD cells in each cluster
table(merged_object_vim_pos_subset$seurat_clusters, merged_object_vim_pos_subset$orig.ident)
#         2020_E15_Rep1   2025_E15.5_Rep1  2025_E15.5_Rep2    COLD_NSC  CON_NSC
# 0             2                241             314            1219    1071
# 1            29               555             621             754     699
# 2            22               237             270             915     419
# 3             5                37              51             198     271
# 4             0                 52            58              41      65
# 5             0                14              16             0      25

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
# 2954     3394    (Total 6348)

saveRDS(clean_combined_object, file = "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/blackshaw + kb/clean_combined_object_blackshaw_and_kb.rds")

# identifying Cycling NSPC cells ------------------------------------------

# cell cycle analysis identified cluster 3 as cycling (done on 'cell cycle workflow on blackshaw+kb cells.R' workflow)

# Cycling NSPCs cell markers
plot <- FeaturePlot(clean_combined_object, features = c("Birc5",
                                                        "Ube2c",
                                                        "Ccnb1",
                                                        "Cdca8",
                                                        "Tpx2"))
ggsave(filename = '/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/blackshaw + kb/analysis on Vim+ subset/Cycling NSPCs.png', plot, width = 15, height = 15)

# cluster 2 is cycling and has cycling NSPC markers 

# naming cluster 2  -------------------------------------------------------
clean_combined_object$cluster_names <- ifelse(
  as.character(clean_combined_object$seurat_clusters) == "2", 
  "Cycling NSPCs", 
  as.character(clean_combined_object$seurat_clusters))

# colour
unique_names <- unique(clean_combined_object$cluster_names)
cluster_colors1 <- scales::hue_pal()(length(unique_names))
names(cluster_colors1) <- unique_names

# Overwrite the "Cycling NSPCs" slot with your specific green color
cluster_colors1["Cycling NSPCs"] <- "#00BA38"

# store in metadata
# If the cell name has "COLD", label it "COLD". Otherwise, label it "CON".
clean_combined_object$treatment <- ifelse(grepl("COLD", colnames(clean_combined_object)), "COLD", "CON")

# If the cell name has "blackshaw", label it "Blackshaw". Otherwise, label it "KB".
clean_combined_object$dataset <- ifelse(grepl("blackshaw", colnames(clean_combined_object)), "Blackshaw", "KB")

# (This assumes orig.ident accurately holds your 5 original sample names)
clean_combined_object$replicate_group <- clean_combined_object$orig.ident

# Plot UMAP using your custom names and colors
DimPlot(
  clean_combined_object, 
  reduction = "umap", 
  group.by = "cluster_names", 
  cols = cluster_colors1
) + ggtitle("Merged Blackshaw + KB NSC Dataset")

table(clean_combined_object$sex)
# Female   Male 
#  2954   3394 

table(clean_combined_object$treatment)
#  COLD  CON 
#  2397 3951 

table(clean_combined_object$replicate_group)
#  2020_E15_Rep1  2025_E15.5_Rep1    2025_E15.5_Rep2        COLD_NSC         CON_NSC 
#  29             823                1000                     2397            2099 

table(clean_combined_object$cluster_names)
#         0             1             3             4             5       Cycling NSPCs 
#      2288          1986           432           164            29          1449 

saveRDS(clean_combined_object, file = "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/blackshaw + kb/clean_combined_object_blackshaw_and_kb_classified_cycling_nspcs.rds")

# load in clean combined object -------------------------------------------

clean_combined_object <- readRDS('/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/blackshaw + kb/clean_combined_object_blackshaw_and_kb_classified_cycling_nspcs.rds')

# This creates a 2-way table showing the male/female split for ALL clusters
table(clean_combined_object$cluster_names, clean_combined_object$sex)

#                 Female  Male
# 0                1086  1202
# 1                856   1130
# 2                207   225
# 4                 70   94
# 5                 13   16
# Cycling NSPCs    722  727


# umaps for thesis --------------------------------------------------------

# all cells umap in clusters
clean_combined_object <- readRDS('/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/blackshaw + kb/clean_combined_object_blackshaw_and_kb.rds')

DimPlot(clean_combined_object, reduction = "umap")

# just cycling nspcs cells 
clean_combined_object <- readRDS('/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/blackshaw + kb/clean_combined_object_blackshaw_and_kb_classified_cycling_nspcs.rds')

# 1. Find the cell barcodes for the Cycling NSPCs using your metadata column
cycling_cells <- WhichCells(clean_combined_object, expression = cluster_names == "Cycling NSPCs")

# 2. Generate the UMAP with the cycling cells highlighted
DimPlot(clean_combined_object, 
        reduction = "umap", 
        cells.highlight = cycling_cells, 
        cols.highlight = "#00BA38", # Your custom green color
        cols = "grey60",            # Background color for all other cells
        pt.size = 0.8) +
  labs(title = "Cycling NSPCs Highlighted") +
  theme(plot.title = element_text(hjust = 0.5))


# cells split by datasets
clean_combined_object <- readRDS('/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/blackshaw + kb/clean_combined_object_blackshaw_and_kb_classified_cycling_nspcs.rds')

# 1. Ensure your highlight group metadata is set up individually for all 4 groups
clean_combined_object$highlight_group <- "Other Cells"
clean_combined_object$highlight_group[clean_combined_object$replicate_group == "CON_NSC"]      <- "KB CON NSC"
clean_combined_object$highlight_group[clean_combined_object$replicate_group == "COLD_NSC"]     <- "KB COLD NSC"
clean_combined_object$highlight_group[clean_combined_object$replicate_group == "2020_E15_Rep1"]   <- "2020 E15 Rep1"
clean_combined_object$highlight_group[clean_combined_object$replicate_group == "2025_E15.5_Rep1"] <- "2025 E15.5 Rep1"
clean_combined_object$highlight_group[clean_combined_object$replicate_group == "2025_E15.5_Rep2"] <- "2025 E15.5 Rep2"

# 2. Define the cell barcodes for each separate category
con_cells      <- WhichCells(clean_combined_object, expression = highlight_group == "KB CON NSC")
cold_cells     <- WhichCells(clean_combined_object, expression = highlight_group == "KB COLD NSC")
e15_rep1_cells <- WhichCells(clean_combined_object, expression = highlight_group == "2020 E15 Rep1")
e5_rep1_cells  <- WhichCells(clean_combined_object, expression = highlight_group == "2025 E15.5 Rep1")
e5_rep2_cells  <- WhichCells(clean_combined_object, expression = highlight_group == "2025 E15.5 Rep2")

# 3. Create individual UMAPs for each distinct group

# --- Plot 1: KB CON NSC (Pink shade 1) ---
DimPlot(clean_combined_object, reduction = "umap", 
              cells.highlight = con_cells, 
              cols.highlight = "#EE6AA7", cols = "grey80", pt.size = 0.5) +
  labs(title = "This Study (CON)") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")

# --- Plot 2: KB COLD NSC (Pink shade 2 / or same pink) ---
DimPlot(clean_combined_object, reduction = "umap", 
              cells.highlight = cold_cells, 
              cols.highlight = "#FF1493", cols = "grey80", pt.size = 0.5) +
  labs(title = "This Study (COLD)") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")

# --- Plot 3: 2020 E15 Rep1 (Purple) ---
DimPlot(clean_combined_object, reduction = "umap", 
              cells.highlight = e15_rep1_cells, 
              cols.highlight = "#9A32CD", cols = "grey80", pt.size = 0.5) +
  labs(title = "E15 Dataset (Replicate 1) Kim et al., 2020") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")

# --- Plot 4: 2025 E15.5 Rep1 (Blue) ---
DimPlot(clean_combined_object, reduction = "umap", 
              cells.highlight = e5_rep1_cells, 
              cols.highlight = "#377EB8", cols = "grey80", pt.size = 0.5) +
  labs(title = "E15.5 Dataset (Replicate 1) Kim et al., 2025") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")

# --- Plot 5: 2025 E15.5 Rep2 (Turquoise) ---
DimPlot(clean_combined_object, reduction = "umap", 
              cells.highlight = e5_rep2_cells, 
              cols.highlight = "#00C5CD", cols = "grey80", pt.size = 0.5) +
  labs(title = "E15.5 Dataset (Replicate 2) Kim et al., 2025") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")



# cells split by sex and treatment
clean_combined_object <- readRDS('/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/blackshaw + kb/clean_combined_object_blackshaw_and_kb.rds')

# 1. Ensure metadata is set up on your main object
clean_combined_object$treatment <- ifelse(grepl("COLD", colnames(clean_combined_object)), "Cold", "Control")
clean_combined_object$sex_treatment <- paste0(clean_combined_object$sex, " - ", clean_combined_object$treatment)

# 2. Subset your main object into the 4 individual conditions
female_con_obj  <- subset(clean_combined_object, sex_treatment == "Female - Control")
female_cold_obj <- subset(clean_combined_object, sex_treatment == "Female - Cold")
male_con_obj    <- subset(clean_combined_object, sex_treatment == "Male - Control")
male_cold_obj   <- subset(clean_combined_object, sex_treatment == "Male - Cold")

# 3. Generate individual UMAPs for each subset object
umap_female_con <- DimPlot(female_con_obj, reduction = "umap", pt.size = 0.5) + 
  ggtitle("Female - Control") + 
  theme(plot.title = element_text(hjust = 0.5))

umap_female_cold <- DimPlot(female_cold_obj, reduction = "umap", pt.size = 0.5) + 
  ggtitle("Female - Cold") + 
  theme(plot.title = element_text(hjust = 0.5))

umap_male_con <- DimPlot(male_con_obj, reduction = "umap", pt.size = 0.5) + 
  ggtitle("Male - Control") + 
  theme(plot.title = element_text(hjust = 0.5))

umap_male_cold <- DimPlot(male_cold_obj, reduction = "umap", pt.size = 0.5) + 
  ggtitle("Male - Cold") + 
  theme(plot.title = element_text(hjust = 0.5))

# 4. Combine them into your multi-panel layout 
(umap_female_con | umap_female_cold | umap_male_con | umap_male_cold)
