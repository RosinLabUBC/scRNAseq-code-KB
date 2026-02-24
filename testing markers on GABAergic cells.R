
library(Seurat)
library(dplyr)
library(patchwork)
library(clustree)
library(pheatmap)
library(EnhancedVolcano)

# load in object
merged_object_gad_pos <- readRDS('/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Gad1+ : Gad2+ cells/objects/merged_object_gad_pos.rds')

# testing parameters on 'testing markers on GABAergic cells.R' workflow
# 35dims, 0.1 res did not look too good

# clustree ----------------------------------------------------------------
# 30 dims
merged_object_gad_pos <- FindNeighbors(merged_object_gad_pos, dims = 1:30)
merged_object_gad_pos <- FindClusters(merged_object_gad_pos, resolution = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.7, 0.9, 1.3))

plot <- clustree(merged_object_gad_pos, prefix = "SCT_snn_res.")
ggsave(filename = '/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Gad1+ : Gad2+ cells/testing parameters on gad+ subset/30_dims.png', plot, width = 10, height = 9)

# umap with different res

merged_object_gad_pos <- RunUMAP(merged_object_gad_pos, dims = 1:30)

plot <- DimPlot(merged_object_gad_pos, reduction = "umap", group.by = "SCT_snn_res.0.1")
ggsave(filename = '/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Gad1+ : Gad2+ cells/testing parameters on gad+ subset/30 dims/30_dims_umap_0.1.png', plot, width = 10, height = 7)

plot <- DimPlot(merged_object_gad_pos, reduction = "umap", group.by = "SCT_snn_res.0.15")
ggsave(filename = '/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Gad1+ : Gad2+ cells/testing parameters on gad+ subset/30 dims/30_dims_umap_0.15.png', plot, width = 10, height = 7)


# 33 dims
merged_object_gad_pos <- FindNeighbors(merged_object_gad_pos, dims = 1:33)
merged_object_gad_pos <- FindClusters(merged_object_gad_pos, resolution = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.7, 0.9, 1.3))

plot <- clustree(merged_object_gad_pos, prefix = "SCT_snn_res.")
ggsave(filename = '/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Gad1+ : Gad2+ cells/testing parameters on gad+ subset/33_dims.png', plot, width = 10, height = 9)

# umap with different res

merged_object_gad_pos <- RunUMAP(merged_object_gad_pos, dims = 1:33)

plot <- DimPlot(merged_object_gad_pos, reduction = "umap", group.by = "SCT_snn_res.0.1")
ggsave(filename = '/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Gad1+ : Gad2+ cells/testing parameters on gad+ subset/33 dims/33_dims_umap_0.1.png', plot, width = 10, height = 7)

plot <- DimPlot(merged_object_gad_pos, reduction = "umap", group.by = "SCT_snn_res.0.15")
ggsave(filename = '/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Gad1+ : Gad2+ cells/testing parameters on gad+ subset/33 dims/33_dims_umap_0.15.png', plot, width = 10, height = 7)


# 35 dims
merged_object_gad_pos <- FindNeighbors(merged_object_gad_pos, dims = 1:35)
merged_object_gad_pos <- FindClusters(merged_object_gad_pos, resolution = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.7, 0.9, 1.3))

plot <- clustree(merged_object_gad_pos, prefix = "SCT_snn_res.")
ggsave(filename = '/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Gad1+ : Gad2+ cells/testing parameters on gad+ subset/35_dims.png', plot, width = 10, height = 9)

# umap with different res

merged_object_gad_pos <- RunUMAP(merged_object_gad_pos, dims = 1:35)

plot <- DimPlot(merged_object_gad_pos, reduction = "umap", group.by = "SCT_snn_res.0.1")
ggsave(filename = '/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Gad1+ : Gad2+ cells/testing parameters on gad+ subset/35 dims/35_dims_umap_0.1.png', plot, width = 10, height = 7)

plot <- DimPlot(merged_object_gad_pos, reduction = "umap", group.by = "SCT_snn_res.0.15")
ggsave(filename = '/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Gad1+ : Gad2+ cells/testing parameters on gad+ subset/35 dims/35_dims_umap_0.15.png', plot, width = 10, height = 7)

plot <- DimPlot(merged_object_gad_pos, reduction = "umap", group.by = "SCT_snn_res.0.2")
ggsave(filename = '/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Gad1+ : Gad2+ cells/testing parameters on gad+ subset/35 dims/35_dims_umap_0.2.png', plot, width = 10, height = 7)


# 37 dims
merged_object_gad_pos <- FindNeighbors(merged_object_gad_pos, dims = 1:37)
merged_object_gad_pos <- FindClusters(merged_object_gad_pos, resolution = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.7, 0.9, 1.3))

plot <- clustree(merged_object_gad_pos, prefix = "SCT_snn_res.")
ggsave(filename = '/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Gad1+ : Gad2+ cells/testing parameters on gad+ subset/37_dims.png', plot, width = 10, height = 9)

# umap with different res

merged_object_gad_pos <- RunUMAP(merged_object_gad_pos, dims = 1:37)

plot <- DimPlot(merged_object_gad_pos, reduction = "umap", group.by = "SCT_snn_res.0.1")
ggsave(filename = '/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Gad1+ : Gad2+ cells/testing parameters on gad+ subset/37 dims/37_dims_umap_0.1.png', plot, width = 10, height = 7)

plot <- DimPlot(merged_object_gad_pos, reduction = "umap", group.by = "SCT_snn_res.0.15")
ggsave(filename = '/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Gad1+ : Gad2+ cells/testing parameters on gad+ subset/37 dims/37_dims_umap_0.15.png', plot, width = 10, height = 7)

plot <- DimPlot(merged_object_gad_pos, reduction = "umap", group.by = "SCT_snn_res.0.2")
ggsave(filename = '/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Gad1+ : Gad2+ cells/testing parameters on gad+ subset/37 dims/37_dims_umap_0.2.png', plot, width = 10, height = 7)


# 40 dims
merged_object_gad_pos <- FindNeighbors(merged_object_gad_pos, dims = 1:40)
merged_object_gad_pos <- FindClusters(merged_object_gad_pos, resolution = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.7, 0.9, 1.3))

plot <- clustree(merged_object_gad_pos, prefix = "SCT_snn_res.")
ggsave(filename = '/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Gad1+ : Gad2+ cells/testing parameters on gad+ subset/40_dims.png', plot, width = 10, height = 9)

# umap with different res

merged_object_gad_pos <- RunUMAP(merged_object_gad_pos, dims = 1:40)

plot <- DimPlot(merged_object_gad_pos, reduction = "umap", group.by = "SCT_snn_res.0.1")
ggsave(filename = '/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Gad1+ : Gad2+ cells/testing parameters on gad+ subset/40 dims/40_dims_umap_0.1.png', plot, width = 10, height = 7)

plot <- DimPlot(merged_object_gad_pos, reduction = "umap", group.by = "SCT_snn_res.0.15")
ggsave(filename = '/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Gad1+ : Gad2+ cells/testing parameters on gad+ subset/40 dims/40_dims_umap_0.15.png', plot, width = 10, height = 7)


# chosen parameters -------------------------------------------------------

merged_object_gad_pos <- FindNeighbors(merged_object_gad_pos, dims = 1:35)	
merged_object_gad_pos <- FindClusters(merged_object_gad_pos, resolution = 0.1)		
merged_object_gad_pos <- RunUMAP(merged_object_gad_pos, dims = 1:35)
DimPlot(merged_object_gad_pos, reduction = "umap")

# testing markers from JR, want to cluster these cells 
# Gad1, Gad2, Sst, Meis2, Dlx1/2/5, Lhx6
# right now using 35dims, 0.1res (same as "merged_object_vim_pos" parameters)

plot <- FeaturePlot(merged_object_gad_pos, features = c("Gad1",
                                                        "Gad2",
                                                        "Sst",
                                                        "Meis2",
                                                        "Dlx1",
                                                        "Dlx2",
                                                        "Dlx5",
                                                        "Lhx6"))
ggsave(filename = '/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Gad1+ : Gad2+ cells/testing parameters on gad+ subset/41 dims/41_dims_0.1_res_feature.png', plot, width = 12, height = 9)


plot <- VlnPlot(merged_object_gad_pos, features = c("Gad1",
                                                        "Gad2",
                                                        "Sst",
                                                        "Meis2",
                                                        "Dlx1",
                                                        "Dlx2",
                                                        "Dlx5",
                                                        "Lhx6"))
ggsave(filename = '/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Gad1+ : Gad2+ cells/testing parameters on gad+ subset/41 dims/41_dims_0.1_res_vln.png', plot, width = 10, height = 10)


# JR: you can get rid of Pvalb and Lamp5 for future analysis (doesn't look like it lights anything up)
# JR: but to try to get better clustering of Sst and Meis2, I will get you to try slightly fewer and slightly more clusters (change dims and res), but as you said, this can be done later

# JR said to go with 35 dims, 0.1 res