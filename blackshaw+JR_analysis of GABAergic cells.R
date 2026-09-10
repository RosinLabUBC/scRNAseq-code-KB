
library(Seurat)
library(dplyr)
library(patchwork)
library(clustree)
library(pheatmap)
library(EnhancedVolcano)

# load in merged_object
merged_object <- readRDS('/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/blackshaw + kb/merged_object_SCTransformed.rds')

# Confirm number of cells
ncol(merged_object)
# 34328 cells

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

FeaturePlot(merged_object, features = c ("Gad1", "Gad2"))

# Subset cells expressing Gad1 or Gad2 or both
merged_object_gad_pos <- subset(merged_object, subset = Gad1 > 0 | Gad2 > 0)

# Confirm number of cells
ncol(merged_object_gad_pos)
# 16351 cells

# Number of cells by condition
table(merged_object_gad_pos$orig.ident)
#    2020_E15_Rep1   2025_E15.5_Rep1   2025_E15.5_Rep2        COLD_NSC         CON_NSC 
#        36              1856            2203                   6853            5403 

# check
FeaturePlot(merged_object_gad_pos, features = c ("Gad1", "Gad2"))

DefaultAssay(merged_object_gad_pos) <- "SCT"
merged_object_gad_pos <- RunPCA(merged_object_gad_pos)
ElbowPlot(merged_object_gad_pos, ndims = 50, reduction = "pca")		

saveRDS(merged_object_gad_pos, "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/blackshaw + kb/merged_object_gad_pos.rds")

# load in Gad+ subsetted cells --------------------------------------------
# tested parameters -- chose the best clustree
# going with 35 dims, 0.1 res --- clustree looked best
# run DEG analysis on all cells 

merged_object_gad_pos <- readRDS('/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/blackshaw + kb/merged_object_gad_pos.rds')

# Confirm number of cells
ncol(merged_object_gad_pos)
# 16351 cells

# going with 35 dims, 0.1 res 
merged_object_gad_pos <- FindNeighbors(merged_object_gad_pos, dims = 1:35)	
merged_object_gad_pos <- FindClusters(merged_object_gad_pos, resolution = 0.1)		
merged_object_gad_pos <- RunUMAP(merged_object_gad_pos, dims = 1:35)
DimPlot(merged_object_gad_pos, reduction = "umap")
DimPlot(merged_object_gad_pos, reduction = "umap", group.by = c("orig.ident", "seurat_clusters"), split.by = "orig.ident")

# removing male and female double positive and non-expressing cells -------
# define markers and expression threshold
female_genes <- c("Xist")
male_genes <- c("Ddx3y", "Uty", "Eif2s3y")
expr_threshold <- 0.1

# Extract SCTransformed data matrix
expr <- GetAssayData(merged_object_gad_pos, assay = "SCT", layer = "data")

# Identify cells expressing male or female markers
female_expr <- colSums(expr[female_genes, , drop = FALSE] > expr_threshold) > 0
male_expr   <- colSums(expr[male_genes, , drop = FALSE] > expr_threshold) > 0

# Assign sex labels: only keep non-overlapping classifications
sex_labels <- ifelse(female_expr & !male_expr, "Female",
                     ifelse(male_expr & !female_expr, "Male", NA))

# Store in metadata
merged_object_gad_pos$sex <- sex_labels

# Directly use the metadata to filter valid cells
valid_cells <- rownames(merged_object_gad_pos@meta.data)[!is.na(merged_object_gad_pos@meta.data$sex)]

# Subset by cell names
merged_object_gad_pos_clean <- subset(merged_object_gad_pos, cells = valid_cells)

# Optional check
table(merged_object_gad_pos_clean$sex)
# Female   Male 
# 6161     6392    (Total 12,553)

# store in metadata
# If the cell name has "COLD", label it "COLD". Otherwise, label it "CON".
merged_object_gad_pos_clean$treatment <- ifelse(grepl("COLD", colnames(merged_object_gad_pos_clean)), "COLD", "CON")

# If the cell name has "blackshaw", label it "Blackshaw". Otherwise, label it "KB".
merged_object_gad_pos_clean$dataset <- ifelse(grepl("blackshaw", colnames(merged_object_gad_pos_clean)), "Blackshaw", "KB")

# (This assumes orig.ident accurately holds your 5 original sample names)
merged_object_gad_pos_clean$replicate_group <- merged_object_gad_pos_clean$orig.ident

saveRDS(merged_object_gad_pos_clean, file = "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/blackshaw + kb/merged_object_gad_pos_clean.rds")
