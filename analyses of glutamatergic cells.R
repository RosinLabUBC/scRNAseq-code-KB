# on all cells (not just Vim+), looking for cells positive for glutamatergic markers 

library(Seurat)
library(dplyr)
library(patchwork)
library(clustree)
library(pheatmap)
library(EnhancedVolcano)

# clustering merged object that is SCTransformed --------------------------
# merged object contains data sets COLD M & F and CON M & F, they have been dieted then SCTransformed again

# load in merged_object
merged_object <- readRDS('/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/merged_object_SCTransformed.rds')

# Confirm number of cells
ncol(merged_object)
# 22643 cells

# Number of cells by condition
table(merged_object$orig.ident)
# COLD_NSC  CON_NSC 
# 12564    10079 

# umap showing cells
DefaultAssay(merged_object) <- "SCT"
merged_object <- RunPCA(merged_object)
ElbowPlot(merged_object, ndims = 50, reduction = "pca")

merged_object <- FindNeighbors(merged_object, dims = 1:35)	
merged_object <- FindClusters(merged_object, resolution = 0.15)		
merged_object <- RunUMAP(merged_object, dims = 1:35)
DimPlot(merged_object, reduction = "umap")
DimPlot(merged_object, reduction = "umap", group.by = c("orig.ident", "seurat_clusters"), split.by = "orig.ident")
DimPlot(merged_object, reduction = "umap", group.by = "orig.ident", label = F)

FeaturePlot(merged_object, features = c ("Slc17a7", "Slc17a6")) # VGlut1 and VGlut2

FeaturePlot(merged_object, features = c ("Grin1", "Grin2b", "Grin2d")) # NMDA subunits

FeaturePlot(merged_object, features = c ("Gria1", "Gria2", "Gria3", "Gria4")) # AMPA subunits

FeaturePlot(merged_object, features = c ("Grik1", "Grik2", "Grik3")) # Kainate subunits

# going with VGlut1 and VGlut2 markers JR "using only vGlu1 and vGlut2 would pull out the cells that produce glutamate"

# Subset cells expressing VGlut1 or VGlut2 or both
merged_object_vglut_pos <- subset(merged_object, subset = Slc17a7 > 0 | Slc17a6 > 0)

# Confirm number of cells
ncol(merged_object_vglut_pos)
# 5944 cells

# Number of cells by condition
table(merged_object_vglut_pos$orig.ident)
# COLD_NSC  CON_NSC 
# 3119    2825 

# check
FeaturePlot(merged_object_vglut_pos, features = c ("Slc17a7", "Slc17a6"))

DefaultAssay(merged_object_vglut_pos) <- "SCT"
merged_object_vglut_pos <- RunPCA(merged_object_vglut_pos)
ElbowPlot(merged_object_vglut_pos, ndims = 50, reduction = "pca")		

saveRDS(merged_object_vglut_pos, "/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis of Glutamatergic Cells/objects/merged_object_vglut_pos.rds")

# load in VGlut+ subsetted cells --------------------------------------------

merged_object_vglut_pos <- readRDS('/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis of Glutamatergic Cells/objects/merged_object_vglut_pos.rds')

# Confirm number of cells
ncol(merged_object_vglut_pos)
# 5944 cells

# testing parameters on 'testing markers on glutamatergic cells.R' workflow
# going with 35 dims, 0.05 res -- confirmed with JR
merged_object_vglut_pos <- FindNeighbors(merged_object_vglut_pos, dims = 1:35)	
merged_object_vglut_pos <- FindClusters(merged_object_vglut_pos, resolution = 0.05)		
merged_object_vglut_pos <- RunUMAP(merged_object_vglut_pos, dims = 1:35)
DimPlot(merged_object_vglut_pos, reduction = "umap")
DimPlot(merged_object_vglut_pos, reduction = "umap", group.by = c("orig.ident", "seurat_clusters"), split.by = "orig.ident")

# removing male and female double positive and non-expressing cells -------
# define markers and expression threshold
female_genes <- c("Xist")
male_genes <- c("Ddx3y", "Uty", "Eif2s3y")
expr_threshold <- 0.1

# Extract SCTransformed data matrix
expr <- GetAssayData(merged_object_vglut_pos, assay = "SCT", layer = "data")

# Identify cells expressing male or female markers
female_expr <- colSums(expr[female_genes, , drop = FALSE] > expr_threshold) > 0
male_expr   <- colSums(expr[male_genes, , drop = FALSE] > expr_threshold) > 0

# Assign sex labels: only keep non-overlapping classifications
sex_labels <- ifelse(female_expr & !male_expr, "Female",
                     ifelse(male_expr & !female_expr, "Male", NA))

# Store in metadata
merged_object_vglut_pos$sex <- sex_labels

# Directly use the metadata to filter valid cells
valid_cells <- rownames(merged_object_vglut_pos@meta.data)[!is.na(merged_object_vglut_pos@meta.data$sex)]

# Subset by cell names
merged_object_vglut_pos_clean <- subset(merged_object_vglut_pos, cells = valid_cells)

# Optional check
table(merged_object_vglut_pos_clean$sex)
# Female   Male 
# 2416    2250    (Total 4666)

# store treatment in metadata
merged_object_vglut_pos_clean$treatment <- ifelse(grepl("COLD", merged_object_vglut_pos_clean$orig.ident), "COLD", "CON")

table(merged_object_vglut_pos_clean$sex)
# Female   Male 
# 2416    2250 

table(merged_object_vglut_pos_clean$treatment)
# COLD  CON 
# 2295 2371 

saveRDS(merged_object_vglut_pos_clean, file = "/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis of Glutamatergic Cells/objects/merged_object_vglut_pos_clean.rds")

# umaps for vglut+ clean cells ----------------------------------------------

merged_object_vglut_pos_clean <- readRDS('/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis of Glutamatergic Cells/objects/merged_object_vglut_pos_clean.rds')

DimPlot(merged_object_vglut_pos_clean, reduction = "umap")
DimPlot(merged_object_vglut_pos_clean, reduction = "umap", group.by = c("orig.ident", "seurat_clusters"), split.by = "orig.ident")

# add colour identities
cluster_colors <- c(
  "0" = "#F8766D",
  "1" = "#a3a500",
  "2" = "#02bf7d",
  "3" = "#01b0f6",
  "4" = "#e76bf3")

# Female + CON
female_con_obj <- subset(merged_object_vglut_pos_clean, subset = sex == "Female" & treatment == "CON")

# Female + COLD
female_cold_obj <- subset(merged_object_vglut_pos_clean, subset = sex == "Female" & treatment == "COLD")

# Male + CON
male_con_obj <- subset(merged_object_vglut_pos_clean, subset = sex == "Male" & treatment == "CON")

# Male + COLD
male_cold_obj <- subset(merged_object_vglut_pos_clean, subset = sex == "Male" & treatment == "COLD")

# Generate UMAPs using DimPlot
umap_female_con <- DimPlot(
  female_con_obj, reduction = "umap", group.by = "seurat_clusters") +
  ggtitle("Female - Control") +
  scale_color_manual(values = cluster_colors)

umap_female_cold <- DimPlot(
  female_cold_obj, reduction = "umap", group.by = "seurat_clusters") +
  ggtitle("Female - Cold") +
  scale_color_manual(values = cluster_colors)

umap_male_con <- DimPlot(
  male_con_obj, reduction = "umap", group.by = "seurat_clusters") +
  ggtitle("Male - Control") +
  scale_color_manual(values = cluster_colors)

umap_male_cold <- DimPlot(
  male_cold_obj, reduction = "umap", group.by = "seurat_clusters") +
  ggtitle("Male - Cold") +
  scale_color_manual(values = cluster_colors)

# Combine UMAPs into a single panel
umap_female_con | umap_female_cold | umap_male_con | umap_male_cold


# get unique gene expression for each cluster -----------------------------

merged_object_vglut_pos_clean <- readRDS('/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis of Glutamatergic Cells/objects/merged_object_vglut_pos_clean.rds')

# Set default assay to RNA
DefaultAssay(merged_object_vglut_pos_clean) <- "RNA"

# Normalize and scale RNA assay
merged_object_vglut_pos_clean <- NormalizeData(merged_object_vglut_pos_clean)
merged_object_vglut_pos_clean <- FindVariableFeatures(merged_object_vglut_pos_clean)
merged_object_vglut_pos_clean <- ScaleData(merged_object_vglut_pos_clean, features = row.names(merged_object_vglut_pos_clean))

# Join Layers
merged_object_vglut_pos_clean <- JoinLayers(merged_object_vglut_pos_clean, layers = "data")

# Set Idents to clusters
Idents(merged_object_vglut_pos_clean) <- "seurat_clusters"

markers <- FindAllMarkers(
  merged_object_vglut_pos_clean,
  assay = "RNA",
  slot = "data",
  only.pos = TRUE,
  min.pct = 0.1,
  logfc.threshold = 0.25)

top30_per_cluster <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 30)

write.csv(
  top30_per_cluster,
  "/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis of Glutamatergic Cells/Analysis on VGlut+ clean cells/unique gene expression of each cluster/top30_markers_per_cluster.csv",
  row.names = FALSE)

# run DEGs for all of the cells -----------------------------------------
# CON M vs CON F and CON F vs COLD F and CON M vs COLD M

# load in object
merged_object_vglut_pos_clean <- readRDS('/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis of Glutamatergic Cells/objects/merged_object_vglut_pos_clean.rds')

# Set default assay to RNA
DefaultAssay(merged_object_vglut_pos_clean) <- "RNA"

# Normalize and scale RNA assay
merged_object_vglut_pos_clean <- NormalizeData(merged_object_vglut_pos_clean)
merged_object_vglut_pos_clean <- FindVariableFeatures(merged_object_vglut_pos_clean)
merged_object_vglut_pos_clean <- ScaleData(merged_object_vglut_pos_clean, features = row.names(merged_object_vglut_pos_clean))

# Join Layers
merged_object_vglut_pos_clean <- JoinLayers(merged_object_vglut_pos_clean, layers = "data")

# Subset to just the CON cells
con_obj <- subset(merged_object_vglut_pos_clean, subset = treatment == "CON")

# Confirm identities
Idents(con_obj) <- "sex"

# DEGs on all CON M vs CON F ----------------------------------------------

# calcualte DEGs
all_CON_DEGs <- FindMarkers(
  con_obj,
  ident.1 = "Male", # comparing how Male genes change compare to female genes
  ident.2 = "Female",
  assay = "RNA",  
  logfc.threshold = 0.25,
  min.pct = 0.1)

# Save
write.csv(all_CON_DEGs, "/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis of Glutamatergic Cells/Analysis on VGlut+ clean cells/DEGs/all_CON_DEGs_Male_vs_Female_(no_p-value_cutoff)_vglut_pos_cells_only.csv")

saveRDS(all_CON_DEGs, "/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis of Glutamatergic Cells/Analysis on VGlut+ clean cells/DEGs/all_CON_DEGs_Male_vs_Female_(no_p-value_cutoff)_vglut_pos_cells_only.rds")

# DEGs on all CON F vs COLD F -----------------------------------------------------

# Subset out female cells
female_obj <- subset(merged_object_vglut_pos_clean, subset = sex == "Female")

# Confirm identities
Idents(female_obj) <- "treatment"

# calcualte DEGs
all_Female_DEGs <- FindMarkers(
  female_obj,
  ident.1 = "COLD", # comparing how Cold genes change compare to CON genes
  ident.2 = "CON",
  assay = "RNA",  
  logfc.threshold = 0.25,
  min.pct = 0.1)

# Save
write.csv(all_Female_DEGs, "/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis of Glutamatergic Cells/Analysis on VGlut+ clean cells/DEGs/all_Female_DEGs_CON_vs_COLD_(no_p-value_cutoff)_vglut_pos_cells_only.csv")

saveRDS(all_Female_DEGs, "/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis of Glutamatergic Cells/Analysis on VGlut+ clean cells/DEGs/all_Female_DEGs_CON_vs_COLD_(no_p-value_cutoff)_vglut_pos_cells_only.rds")


# DEGs on all CON M vs COLD M -----------------------------------------------------

# Subset out female cells
male_obj <- subset(merged_object_vglut_pos_clean, subset = sex == "Male")

# Confirm identities
Idents(male_obj) <- "treatment"

# calcualte DEGs
all_Male_DEGs <- FindMarkers(
  male_obj,
  ident.1 = "COLD", # comparing how Cold genes change compare to CON genes
  ident.2 = "CON",
  assay = "RNA",  
  logfc.threshold = 0.25,
  min.pct = 0.1)

# Save
write.csv(all_Male_DEGs, "/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis of Glutamatergic Cells/Analysis on VGlut+ clean cells/DEGs/all_Male_DEGs_CON_vs_COLD_(no_p-value_cutoff)_vglut_pos_cells_only.csv")

saveRDS(all_Male_DEGs, "/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis of Glutamatergic Cells/Analysis on VGlut+ clean cells/DEGs/all_Male_DEGs_CON_vs_COLD_(no_p-value_cutoff)_vglut_pos_cells_only.rds")


# make volcano plots ------------------------------------------------------

# load in DEG object

filepath <- paste0("/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis of Glutamatergic Cells/Analysis on VGlut+ clean cells/DEGs/all_CON_DEGs_Male_vs_Female_(no_p-value_cutoff)_vglut_pos_cells_only.rds")

DEGs <- readRDS(filepath)

# Add gene names
DEGs$gene <- rownames(DEGs)

# Categorize differential expression
DEGs$diffexpressed <- "NO"
DEGs$diffexpressed[DEGs$avg_log2FC > 0 & DEGs$p_val_adj < 0.05] <- "UP"
DEGs$diffexpressed[DEGs$avg_log2FC < 0 & DEGs$p_val_adj < 0.05] <- "DOWN"

# Plot
ggplot(DEGs, aes(x = avg_log2FC, y = -log10(p_val_adj), color = diffexpressed)) +
  geom_point(size = 2) +
  geom_vline(xintercept = c(0, 0), linetype = "dashed", color = "gray") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray") +
  scale_color_manual(values = c("DOWN" = "#00AFBB", "NO" = "gray", "UP" = "darkred"),
                     labels = c("Downregulated", "Not significant", "Upregulated")) +
  theme_minimal() +
  labs(title = paste("Volcano Plot – CON M vs CON F (VGlut1+/VGlut2+ cells)"),
       x = "Log2 Fold Change",
       y = "-log10(Adjusted p-value)",
       color = "Expression") +
  geom_text_repel(data = subset(DEGs, diffexpressed != "NO"),
                  aes(label = gene), size = 3, max.overlaps = 35) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(10, 60, 10, 10))

# CHECK -------------------------------------------------------------------

# All upregulated
upregulated <- DEGs %>%
  filter(avg_log2FC > 0, p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC))

write.csv(upregulated,
          paste0("/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis of Glutamatergic Cells/Analysis on VGlut+ clean cells/volcano plots/check/all_upregulated_CON_M_vs_F.csv"),
          row.names = FALSE)

# All downregulated
downregulated <- DEGs %>%
  filter(avg_log2FC < 0, p_val_adj < 0.05) %>%
  arrange(avg_log2FC) 

write.csv(downregulated,
          paste0("/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis of Glutamatergic Cells/Analysis on VGlut+ clean cells/volcano plots/check/all_downregulated_CON_M_vs_F.csv"),
          row.names = FALSE)

# CON F vs COLD F ---------------------------------------------------------
# load in DEG object

filepath <- paste0("/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis of Glutamatergic Cells/Analysis on VGlut+ clean cells/DEGs/all_Female_DEGs_CON_vs_COLD_(no_p-value_cutoff)_vglut_pos_cells_only.rds")

DEGs <- readRDS(filepath)

# Add gene names
DEGs$gene <- rownames(DEGs)

# Categorize differential expression
DEGs$diffexpressed <- "NO"
DEGs$diffexpressed[DEGs$avg_log2FC > 0 & DEGs$p_val_adj < 0.05] <- "UP"
DEGs$diffexpressed[DEGs$avg_log2FC < 0 & DEGs$p_val_adj < 0.05] <- "DOWN"

# Plot
ggplot(DEGs, aes(x = avg_log2FC, y = -log10(p_val_adj), color = diffexpressed)) +
  geom_point(size = 2) +
  geom_vline(xintercept = c(0, 0), linetype = "dashed", color = "gray") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray") +
  scale_color_manual(values = c("DOWN" = "#00AFBB", "NO" = "gray", "UP" = "darkred"),
                     labels = c("Downregulated", "Not significant", "Upregulated")) +
  theme_minimal() +
  labs(title = paste("Volcano Plot – COLD F vs CON F (VGlut1+/VGlut2+ cells)"),
       x = "Log2 Fold Change",
       y = "-log10(Adjusted p-value)",
       color = "Expression") +
  geom_text_repel(data = subset(DEGs, diffexpressed != "NO"),
                  aes(label = gene), size = 3, max.overlaps = 20) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(10, 60, 10, 10))


# CHECK -------------------------------------------------------------------

# All upregulated
upregulated <- DEGs %>%
  filter(avg_log2FC > 0, p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC))

write.csv(upregulated,
          paste0("/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis of Glutamatergic Cells/Analysis on VGlut+ clean cells/volcano plots/check/all_upregulated_CON_F_vs_COLD_F.csv"),
          row.names = FALSE)

# All downregulated
downregulated <- DEGs %>%
  filter(avg_log2FC < 0, p_val_adj < 0.05) %>%
  arrange(avg_log2FC) 

write.csv(downregulated,
          paste0("/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis of Glutamatergic Cells/Analysis on VGlut+ clean cells/volcano plots/check/all_downregulated_CON_F_vs_COLD_F.csv"),
          row.names = FALSE)


# CON M vs COLD M ---------------------------------------------------------
# load in DEG object

filepath <- paste0("/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis of Glutamatergic Cells/Analysis on VGlut+ clean cells/DEGs/all_Male_DEGs_CON_vs_COLD_(no_p-value_cutoff)_vglut_pos_cells_only.rds")

DEGs <- readRDS(filepath)

# Add gene names
DEGs$gene <- rownames(DEGs)

# Categorize differential expression
DEGs$diffexpressed <- "NO"
DEGs$diffexpressed[DEGs$avg_log2FC > 0 & DEGs$p_val_adj < 0.05] <- "UP"
DEGs$diffexpressed[DEGs$avg_log2FC < 0 & DEGs$p_val_adj < 0.05] <- "DOWN"

# Plot
ggplot(DEGs, aes(x = avg_log2FC, y = -log10(p_val_adj), color = diffexpressed)) +
  geom_point(size = 2) +
  geom_vline(xintercept = c(0, 0), linetype = "dashed", color = "gray") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray") +
  scale_color_manual(values = c("DOWN" = "#00AFBB", "NO" = "gray", "UP" = "darkred"),
                     labels = c("Downregulated", "Not significant", "Upregulated")) +
  theme_minimal() +
  labs(title = paste("Volcano Plot – COLD M vs CON M (VGlut1+/VGlut2+ cells)"),
       x = "Log2 Fold Change",
       y = "-log10(Adjusted p-value)",
       color = "Expression") +
  geom_text_repel(data = subset(DEGs, diffexpressed != "NO"),
                  aes(label = gene), size = 3, max.overlaps = 20) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(10, 60, 10, 10))


# CHECK -------------------------------------------------------------------

# All upregulated
upregulated <- DEGs %>%
  filter(avg_log2FC > 0, p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC))

write.csv(upregulated,
          paste0("/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis of Glutamatergic Cells/Analysis on VGlut+ clean cells/volcano plots/check/all_upregulated_CON_M_vs_COLD_M.csv"),
          row.names = FALSE)

# All downregulated
downregulated <- DEGs %>%
  filter(avg_log2FC < 0, p_val_adj < 0.05) %>%
  arrange(avg_log2FC) 

write.csv(downregulated,
          paste0("/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis of Glutamatergic Cells/Analysis on VGlut+ clean cells/volcano plots/check/all_downregulated_CON_M_vs_COLD_M.csv"),
          row.names = FALSE)

