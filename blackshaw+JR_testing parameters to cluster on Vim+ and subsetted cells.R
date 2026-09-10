# 30 dims -----------------------------------------------------------------

merged_object_vim_pos_subset <- FindNeighbors(merged_object_vim_pos_subset, dims = 1:30)
merged_object_vim_pos_subset <- FindClusters(merged_object_vim_pos_subset, resolution = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.7, 0.9, 1.3))

plot <- clustree(merged_object_vim_pos_subset, prefix = "SCT_snn_res.")
ggsave(filename = '/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/Vim+ cells/35 dims, 0.05 res/cluster 1, Ascl1+/clustree/30_dims.png', plot, width = 10, height = 9)

# 33 dims -----------------------------------------------------------------

merged_object_vim_pos_subset <- FindNeighbors(merged_object_vim_pos_subset, dims = 1:33)
merged_object_vim_pos_subset <- FindClusters(merged_object_vim_pos_subset, resolution = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.7, 0.9, 1.3))

plot <- clustree(merged_object_vim_pos_subset, prefix = "SCT_snn_res.")
ggsave(filename = '/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/Vim+ cells/35 dims, 0.05 res/cluster 1, Ascl1+/clustree/33_dims.png', plot, width = 10, height = 9)

# 35 dims -----------------------------------------------------------------

merged_object_vim_pos_subset <- FindNeighbors(merged_object_vim_pos_subset, dims = 1:35)
merged_object_vim_pos_subset <- FindClusters(merged_object_vim_pos_subset, resolution = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.7, 0.9, 1.3))

plot <- clustree(merged_object_vim_pos_subset, prefix = "SCT_snn_res.")
ggsave(filename = '/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/Vim+ cells/35 dims, 0.05 res/cluster 1, Ascl1+/clustree/35_dims.png', plot, width = 10, height = 9)


# 37 dims -----------------------------------------------------------------

merged_object_vim_pos_subset <- FindNeighbors(merged_object_vim_pos_subset, dims = 1:37)
merged_object_vim_pos_subset <- FindClusters(merged_object_vim_pos_subset, resolution = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.7, 0.9, 1.3))

plot <- clustree(merged_object_vim_pos_subset, prefix = "SCT_snn_res.")
ggsave(filename = '/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/Vim+ cells/35 dims, 0.05 res/cluster 1, Ascl1+/clustree/37_dims.png', plot, width = 10, height = 9)

# 40 dims -----------------------------------------------------------------

merged_object_vim_pos_subset <- FindNeighbors(merged_object_vim_pos_subset, dims = 1:40)
merged_object_vim_pos_subset <- FindClusters(merged_object_vim_pos_subset, resolution = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.7, 0.9, 1.3))

plot <- clustree(merged_object_vim_pos_subset, prefix = "SCT_snn_res.")
ggsave(filename = '/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/Vim+ cells/35 dims, 0.05 res/cluster 1, Ascl1+/clustree/40_dims.png', plot, width = 10, height = 9)



# using 30 dims, 0.2 res -----------------------------
merged_object_vim_pos_subset <- FindNeighbors(merged_object_vim_pos_subset, dims = 1:30)	
merged_object_vim_pos_subset <- FindClusters(merged_object_vim_pos_subset, resolution = 0.2)		
merged_object_vim_pos_subset <- RunUMAP(merged_object_vim_pos_subset, dims = 1:30)
DimPlot(merged_object_vim_pos_subset, reduction = "umap")
DimPlot(merged_object_vim_pos_subset, reduction = "umap", group.by = c("orig.ident", "seurat_clusters"), split.by = "orig.ident")


# using 33 dims, 0.2 res -----------------------------
merged_object_vim_pos_subset <- FindNeighbors(merged_object_vim_pos_subset, dims = 1:33)	
merged_object_vim_pos_subset <- FindClusters(merged_object_vim_pos_subset, resolution = 0.2)		
merged_object_vim_pos_subset <- RunUMAP(merged_object_vim_pos_subset, dims = 1:33)
DimPlot(merged_object_vim_pos_subset, reduction = "umap")
DimPlot(merged_object_vim_pos_subset, reduction = "umap", group.by = c("orig.ident", "seurat_clusters"), split.by = "orig.ident")


