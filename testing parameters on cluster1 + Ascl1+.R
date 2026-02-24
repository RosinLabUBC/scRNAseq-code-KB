# Testing parameters on Vim+, cluster 1, and Ascl+ subset 

# 30 dims -----------------------------------------------------------------

merged_object_vim_pos_subset <- FindNeighbors(merged_object_vim_pos_subset, dims = 1:30)
merged_object_vim_pos_subset <- FindClusters(merged_object_vim_pos_subset, resolution = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.7, 0.9, 1.3))

plot <- clustree(merged_object_vim_pos_subset, prefix = "SCT_snn_res.")
ggsave(filename = '/Users/kodeebao/Documents/Final Analysis v_3/Strategy 2/Vim+ cells only/35 dims, 0.1 res/cluster 1, Ascl1+/clustree/30_dims.png', plot, width = 10, height = 9)

# umap with different res

merged_object_vim_pos_subset <- RunUMAP(merged_object_vim_pos_subset, dims = 1:30)

plot <- DimPlot(merged_object_vim_pos_subset, reduction = "umap", group.by = "SCT_snn_res.0.2")
ggsave(filename = '/Users/kodeebao/Documents/Final Analysis v_3/Strategy 2/Vim+ cells only/35 dims, 0.1 res/cluster 1, Ascl1+/clustree/30_dims_umap_0.2.png', plot, width = 10, height = 7)


# 33 dims -----------------------------------------------------------------

merged_object_vim_pos_subset <- FindNeighbors(merged_object_vim_pos_subset, dims = 1:33)
merged_object_vim_pos_subset <- FindClusters(merged_object_vim_pos_subset, resolution = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.7, 0.9, 1.3))

plot <- clustree(merged_object_vim_pos_subset, prefix = "SCT_snn_res.")
ggsave(filename = '/Users/kodeebao/Documents/Final Analysis v_3/Strategy 2/Vim+ cells only/35 dims, 0.1 res/cluster 1, Ascl1+/clustree/33_dims.png', plot, width = 10, height = 9)

# 35 dims -----------------------------------------------------------------

merged_object_vim_pos_subset <- FindNeighbors(merged_object_vim_pos_subset, dims = 1:35)
merged_object_vim_pos_subset <- FindClusters(merged_object_vim_pos_subset, resolution = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.7, 0.9, 1.3))

plot <- clustree(merged_object_vim_pos_subset, prefix = "SCT_snn_res.")
ggsave(filename = '/Users/kodeebao/Documents/Final Analysis v_3/Strategy 2/Vim+ cells only/35 dims, 0.1 res/cluster 1, Ascl1+/clustree/35_dims.png', plot, width = 10, height = 9)

# umap with different res

merged_object_vim_pos_subset <- RunUMAP(merged_object_vim_pos_subset, dims = 1:35)

plot <- DimPlot(merged_object_vim_pos_subset, reduction = "umap", group.by = "SCT_snn_res.0.1")
ggsave(filename = '/Users/kodeebao/Documents/Final Analysis v_3/Strategy 2/Vim+ cells only/35 dims, 0.1 res/cluster 1, Ascl1+/clustree/35_dims_umap_0.1.png', plot, width = 10, height = 7)

plot <- DimPlot(merged_object_vim_pos_subset, reduction = "umap", group.by = "SCT_snn_res.0.05")
ggsave(filename = '/Users/kodeebao/Documents/Final Analysis v_3/Strategy 2/Vim+ cells only/35 dims, 0.1 res/cluster 1, Ascl1+/clustree/35_dims_umap_0.05.png', plot, width = 10, height = 7)


# 37 dims -----------------------------------------------------------------

merged_object_vim_pos_subset <- FindNeighbors(merged_object_vim_pos_subset, dims = 1:37)
merged_object_vim_pos_subset <- FindClusters(merged_object_vim_pos_subset, resolution = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.7, 0.9, 1.3))

plot <- clustree(merged_object_vim_pos_subset, prefix = "SCT_snn_res.")
ggsave(filename = '/Users/kodeebao/Documents/Final Analysis v_3/Strategy 2/Vim+ cells only/35 dims, 0.1 res/cluster 1, Ascl1+/clustree/37_dims.png', plot, width = 10, height = 9)

# 40 dims -----------------------------------------------------------------

merged_object_vim_pos_subset <- FindNeighbors(merged_object_vim_pos_subset, dims = 1:40)
merged_object_vim_pos_subset <- FindClusters(merged_object_vim_pos_subset, resolution = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.7, 0.9, 1.3))

plot <- clustree(merged_object_vim_pos_subset, prefix = "SCT_snn_res.")
ggsave(filename = '/Users/kodeebao/Documents/Final Analysis v_3/Strategy 2/Vim+ cells only/35 dims, 0.1 res/cluster 1, Ascl1+/clustree/40_dims.png', plot, width = 10, height = 9)


# testing harmonys markers on 35 dims, 0.1 res (merged and Vim+, cluster 1, Ascl1+ subsetted) ------------------

merged_object_vim_pos_subset <- FindNeighbors(merged_object_vim_pos_subset, dims = 1:35)	
merged_object_vim_pos_subset <- FindClusters(merged_object_vim_pos_subset, resolution = 0.1)		
merged_object_vim_pos_subset <- RunUMAP(merged_object_vim_pos_subset, dims = 1:35)
DimPlot(merged_object_vim_pos_subset, reduction = "umap")
DimPlot(merged_object_vim_pos_subset, reduction = "umap", group.by = c("orig.ident", "seurat_clusters"), split.by = "orig.ident")

# proliferative cell markers --------------------------------------------
plot <- FeaturePlot(merged_object_vim_pos_subset, features = c("Ube2c"))
ggsave(filename = '/Users/kodeebao/Documents/Final Analysis v_3/Strategy 2/Vim+ cells only/35 dims, 0.1 res/cluster 1, Ascl1+/35 dims, 0.1 res/proliferative_feature.png', plot, width = 8, height = 5)

plot <- VlnPlot(merged_object_vim_pos_subset, features = c("Ube2c"))
ggsave(filename = '/Users/kodeebao/Documents/Final Analysis v_3/Strategy 2/Vim+ cells only/35 dims, 0.1 res/cluster 1, Ascl1+/35 dims, 0.1 res/proliferative_vln.png', plot, width = 8, height = 5)


# RGCs-glia cell markers --------------------------------------------
plot <- FeaturePlot(merged_object_vim_pos_subset, features = c("Gdpd2"))
ggsave(filename = '/Users/kodeebao/Documents/Final Analysis v_3/Strategy 2/Vim+ cells only/35 dims, 0.1 res/cluster 1, Ascl1+/35 dims, 0.1 res/RGCs-glia_feature.png', plot, width = 8, height = 5)

plot <- VlnPlot(merged_object_vim_pos_subset, features = c("Gdpd2"))
ggsave(filename = '/Users/kodeebao/Documents/Final Analysis v_3/Strategy 2/Vim+ cells only/35 dims, 0.1 res/cluster 1, Ascl1+/35 dims, 0.1 res/RGCs-glia_vln.png', plot, width = 8, height = 5)


# OPCs cell markers --------------------------------------------
plot <- FeaturePlot(merged_object_vim_pos_subset, features = c("Olig2",
                                                        "Sox10",
                                                        "Pdgfra"))
ggsave(filename = '/Users/kodeebao/Documents/Final Analysis v_3/Strategy 2/Vim+ cells only/35 dims, 0.1 res/cluster 1, Ascl1+/35 dims, 0.1 res/OPCs_feature.png', plot, width = 20, height = 10)

plot <- VlnPlot(merged_object_vim_pos_subset, features = c("Olig2",
                                                    "Sox10",
                                                    "Pdgfra"))
ggsave(filename = '/Users/kodeebao/Documents/Final Analysis v_3/Strategy 2/Vim+ cells only/35 dims, 0.1 res/cluster 1, Ascl1+/35 dims, 0.1 res/OPCs_vln.png', plot, width = 20, height = 10)


# Astrocytes cell markers --------------------------------------------
plot <- FeaturePlot(merged_object_vim_pos_subset, features = c("Aldh1l1",
                                                        "Gfap",
                                                        "Slc1a3"))
ggsave(filename = '/Users/kodeebao/Documents/Final Analysis v_3/Strategy 2/Vim+ cells only/35 dims, 0.1 res/cluster 1, Ascl1+/35 dims, 0.1 res/astrocytes_feature.png', plot, width = 20, height = 10)

plot <- VlnPlot(merged_object_vim_pos_subset, features = c("Aldh1l1",
                                                    "Gfap",
                                                    "Slc1a3"))
ggsave(filename = '/Users/kodeebao/Documents/Final Analysis v_3/Strategy 2/Vim+ cells only/35 dims, 0.1 res/cluster 1, Ascl1+/35 dims, 0.1 res/astrocytes_vln.png', plot, width = 20, height = 10)

# blackshaw NPC markers cell markers --------------------------------------------
plot <- FeaturePlot(merged_object_vim_pos_subset, features = c("Ascl1",
                                                        "Neurog2"))
ggsave(filename = '/Users/kodeebao/Documents/Final Analysis v_3/Strategy 2/Vim+ cells only/35 dims, 0.1 res/cluster 1, Ascl1+/35 dims, 0.1 res/blackshaw_NPCs_feature.png', plot, width = 15, height = 5)

plot <- VlnPlot(merged_object_vim_pos_subset, features = c("Ascl1",
                                                    "Neurog2"))
ggsave(filename = '/Users/kodeebao/Documents/Final Analysis v_3/Strategy 2/Vim+ cells only/35 dims, 0.1 res/cluster 1, Ascl1+/35 dims, 0.1 res/blackshaw_NPCs_vln.png', plot, width = 15, height = 5)


# testing harmonys markers on 30 dims, 0.2 res (merged and Vim+, cluster 1, Ascl1+ subsetted) ------------------

merged_object_vim_pos_subset <- FindNeighbors(merged_object_vim_pos_subset, dims = 1:30)	
merged_object_vim_pos_subset <- FindClusters(merged_object_vim_pos_subset, resolution = 0.2)		
merged_object_vim_pos_subset <- RunUMAP(merged_object_vim_pos_subset, dims = 1:30)
DimPlot(merged_object_vim_pos_subset, reduction = "umap")
DimPlot(merged_object_vim_pos_subset, reduction = "umap", group.by = c("orig.ident", "seurat_clusters"), split.by = "orig.ident")

# proliferative cell markers --------------------------------------------
plot <- FeaturePlot(merged_object_vim_pos_subset, features = c("Ube2c"))
ggsave(filename = '/Users/kodeebao/Documents/Final Analysis v_3/Strategy 2/Vim+ cells only/35 dims, 0.1 res/cluster 1, Ascl1+/30 dims, 0.2 res/proliferative_feature.png', plot, width = 8, height = 5)

plot <- VlnPlot(merged_object_vim_pos_subset, features = c("Ube2c"))
ggsave(filename = '/Users/kodeebao/Documents/Final Analysis v_3/Strategy 2/Vim+ cells only/35 dims, 0.1 res/cluster 1, Ascl1+/30 dims, 0.2 res/proliferative_vln.png', plot, width = 8, height = 5)


# RGCs-glia cell markers --------------------------------------------
plot <- FeaturePlot(merged_object_vim_pos_subset, features = c("Gdpd2"))
ggsave(filename = '/Users/kodeebao/Documents/Final Analysis v_3/Strategy 2/Vim+ cells only/35 dims, 0.1 res/cluster 1, Ascl1+/30 dims, 0.2 res/RGCs-glia_feature.png', plot, width = 8, height = 5)

plot <- VlnPlot(merged_object_vim_pos_subset, features = c("Gdpd2"))
ggsave(filename = '/Users/kodeebao/Documents/Final Analysis v_3/Strategy 2/Vim+ cells only/35 dims, 0.1 res/cluster 1, Ascl1+/30 dims, 0.2 res/RGCs-glia_vln.png', plot, width = 8, height = 5)


# OPCs cell markers --------------------------------------------
plot <- FeaturePlot(merged_object_vim_pos_subset, features = c("Olig2",
                                                               "Sox10",
                                                               "Pdgfra"))
ggsave(filename = '/Users/kodeebao/Documents/Final Analysis v_3/Strategy 2/Vim+ cells only/35 dims, 0.1 res/cluster 1, Ascl1+/30 dims, 0.2 res/OPCs_feature.png', plot, width = 20, height = 10)

plot <- VlnPlot(merged_object_vim_pos_subset, features = c("Olig2",
                                                           "Sox10",
                                                           "Pdgfra"))
ggsave(filename = '/Users/kodeebao/Documents/Final Analysis v_3/Strategy 2/Vim+ cells only/35 dims, 0.1 res/cluster 1, Ascl1+/30 dims, 0.2 res/OPCs_vln.png', plot, width = 20, height = 10)


# Astrocytes cell markers --------------------------------------------
plot <- FeaturePlot(merged_object_vim_pos_subset, features = c("Aldh1l1",
                                                               "Gfap",
                                                               "Slc1a3"))
ggsave(filename = '/Users/kodeebao/Documents/Final Analysis v_3/Strategy 2/Vim+ cells only/35 dims, 0.1 res/cluster 1, Ascl1+/30 dims, 0.2 res/astrocytes_feature.png', plot, width = 20, height = 10)

plot <- VlnPlot(merged_object_vim_pos_subset, features = c("Aldh1l1",
                                                           "Gfap",
                                                           "Slc1a3"))
ggsave(filename = '/Users/kodeebao/Documents/Final Analysis v_3/Strategy 2/Vim+ cells only/35 dims, 0.1 res/cluster 1, Ascl1+/30 dims, 0.2 res/astrocytes_vln.png', plot, width = 20, height = 10)

# blackshaw NPC markers cell markers --------------------------------------------
plot <- FeaturePlot(merged_object_vim_pos_subset, features = c("Ascl1",
                                                               "Neurog2"))
ggsave(filename = '/Users/kodeebao/Documents/Final Analysis v_3/Strategy 2/Vim+ cells only/35 dims, 0.1 res/cluster 1, Ascl1+/30 dims, 0.2 res/blackshaw_NPCs_feature.png', plot, width = 15, height = 5)

plot <- VlnPlot(merged_object_vim_pos_subset, features = c("Ascl1",
                                                           "Neurog2"))
ggsave(filename = '/Users/kodeebao/Documents/Final Analysis v_3/Strategy 2/Vim+ cells only/35 dims, 0.1 res/cluster 1, Ascl1+/30 dims, 0.2 res/blackshaw_NPCs_vln.png', plot, width = 15, height = 5)



# -------------------------------------------------------------------------


