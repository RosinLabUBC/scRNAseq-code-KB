
library(Seurat)
library(dplyr)
library(patchwork)
library(clustree)
library(pheatmap)
library(EnhancedVolcano)

# load in object
merged_object_vglut_pos <- readRDS('/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis of Glutamatergic Cells/objects/merged_object_vglut_pos.rds')

# testing parameters on 'testing markers on glutamatergic cells.R' workflow
# clustree ----------------------------------------------------------------
# 30 dims
merged_object_vglut_pos <- FindNeighbors(merged_object_vglut_pos, dims = 1:30)
merged_object_vglut_pos <- FindClusters(merged_object_vglut_pos, resolution = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.7, 0.9, 1.3))

plot <- clustree(merged_object_vglut_pos, prefix = "SCT_snn_res.")
ggsave(filename = '/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis of Glutamatergic Cells/testing parameters on vglut+ subset/30_dims.png', plot, width = 10, height = 9)

# umap with different res

merged_object_vglut_pos <- RunUMAP(merged_object_vglut_pos, dims = 1:30)

plot <- DimPlot(merged_object_vglut_pos, reduction = "umap", group.by = "SCT_snn_res.0.1")
ggsave(filename = '/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis of Glutamatergic Cells/testing parameters on vglut+ subset/30 dims/30_dims_umap_0.1.png', plot, width = 10, height = 7)

plot <- DimPlot(merged_object_vglut_pos, reduction = "umap", group.by = "SCT_snn_res.0.15")
ggsave(filename = '/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis of Glutamatergic Cells/testing parameters on vglut+ subset/30 dims/30_dims_umap_0.15.png', plot, width = 10, height = 7)

# 33 dims
merged_object_vglut_pos <- FindNeighbors(merged_object_vglut_pos, dims = 1:33)
merged_object_vglut_pos <- FindClusters(merged_object_vglut_pos, resolution = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.7, 0.9, 1.3))

plot <- clustree(merged_object_vglut_pos, prefix = "SCT_snn_res.")
ggsave(filename = '/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis of Glutamatergic Cells/testing parameters on vglut+ subset/33_dims.png', plot, width = 10, height = 9)

# umap with different res

# 35 dims
merged_object_vglut_pos <- FindNeighbors(merged_object_vglut_pos, dims = 1:35)
merged_object_vglut_pos <- FindClusters(merged_object_vglut_pos, resolution = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.7, 0.9, 1.3))

plot <- clustree(merged_object_vglut_pos, prefix = "SCT_snn_res.")
ggsave(filename = '/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis of Glutamatergic Cells/testing parameters on vglut+ subset/35_dims.png', plot, width = 10, height = 9)


# 37 dims
merged_object_vglut_pos <- FindNeighbors(merged_object_vglut_pos, dims = 1:37)
merged_object_vglut_pos <- FindClusters(merged_object_vglut_pos, resolution = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.7, 0.9, 1.3))

plot <- clustree(merged_object_vglut_pos, prefix = "SCT_snn_res.")
ggsave(filename = '/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis of Glutamatergic Cells/testing parameters on vglut+ subset/37_dims.png', plot, width = 10, height = 9)


# 40 dims
merged_object_vglut_pos <- FindNeighbors(merged_object_vglut_pos, dims = 1:40)
merged_object_vglut_pos <- FindClusters(merged_object_vglut_pos, resolution = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.7, 0.9, 1.3))

plot <- clustree(merged_object_vglut_pos, prefix = "SCT_snn_res.")
ggsave(filename = '/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis of Glutamatergic Cells/testing parameters on vglut+ subset/40_dims.png', plot, width = 10, height = 9)
