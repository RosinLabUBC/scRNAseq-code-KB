# Following quality control, CON and COLD Seurat objects were merged into a single dataset. Vim+ cells were then subsetted from the merged object, and an initial round of clustering and testing parameters was performed specifically on the Vim+ population 

# clustering COLD M & F and CON M & F data sets, only Vim+ cells

# 30 dims -----------------------------------------------------------------

merged_object_vim_pos <- FindNeighbors(merged_object_vim_pos, dims = 1:30)
merged_object_vim_pos <- FindClusters(merged_object_vim_pos, resolution = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.7, 0.9, 1.3))

plot <- clustree(merged_object_vim_pos, prefix = "SCT_snn_res.")
ggsave(filename = '/Users/kodeebao/Documents/Final Analysis v_2/COLD and CON cells, Vim+/clustree/30_dims.png', plot, width = 10, height = 9)

# 33 dims -----------------------------------------------------------------

merged_object_vim_pos <- FindNeighbors(merged_object_vim_pos, dims = 1:33)
merged_object_vim_pos <- FindClusters(merged_object_vim_pos, resolution = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.7, 0.9, 1.3))

plot <- clustree(merged_object_vim_pos, prefix = "SCT_snn_res.")
ggsave(filename = '/Users/kodeebao/Documents/Final Analysis v_2/COLD and CON cells, Vim+/clustree/33_dims.png', plot, width = 10, height = 9)

# 35 dims -----------------------------------------------------------------

merged_object_vim_pos <- FindNeighbors(merged_object_vim_pos, dims = 1:35)
merged_object_vim_pos <- FindClusters(merged_object_vim_pos, resolution = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.7, 0.9, 1.3))

plot <- clustree(merged_object_vim_pos, prefix = "SCT_snn_res.")
ggsave(filename = '/Users/kodeebao/Documents/Final Analysis v_2/COLD and CON cells, Vim+/clustree/35_dims.png', plot, width = 10, height = 9)

# umap with different res

merged_object_vim_pos <- RunUMAP(merged_object_vim_pos, dims = 1:35)

plot <- DimPlot(merged_object_vim_pos, reduction = "umap", group.by = "SCT_snn_res.0.05")
ggsave(filename = '/Users/kodeebao/Documents/Final Analysis v_2/COLD and CON cells, Vim+/clustree/35_dims_umap_0.05.png', plot, width = 10, height = 7)

plot <- DimPlot(merged_object_vim_pos, reduction = "umap", group.by = "SCT_snn_res.0.1")
ggsave(filename = '/Users/kodeebao/Documents/Final Analysis v_2/COLD and CON cells, Vim+/clustree/35_dims_umap_0.1.png', plot, width = 10, height = 7)


# 37 dims -----------------------------------------------------------------

merged_object_vim_pos <- FindNeighbors(merged_object_vim_pos, dims = 1:37)
merged_object_vim_pos <- FindClusters(merged_object_vim_pos, resolution = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.7, 0.9, 1.3))

plot <- clustree(merged_object_vim_pos, prefix = "SCT_snn_res.")
ggsave(filename = '/Users/kodeebao/Documents/Final Analysis v_2/COLD and CON cells, Vim+/clustree/37_dims.png', plot, width = 10, height = 9)

# 40 dims -----------------------------------------------------------------

merged_object_vim_pos <- FindNeighbors(merged_object_vim_pos, dims = 1:40)
merged_object_vim_pos <- FindClusters(merged_object_vim_pos, resolution = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.7, 0.9, 1.3))

plot <- clustree(merged_object_vim_pos, prefix = "SCT_snn_res.")
ggsave(filename = '/Users/kodeebao/Documents/Final Analysis v_2/COLD and CON cells, Vim+/clustree/40_dims.png', plot, width = 10, height = 9)

# testing harmonys markers on 35 dims, 0.05 res (merged and Vim+ subsetted) ------------------------

merged_object_vim_pos <- FindNeighbors(merged_object_vim_pos, dims = 1:35)	
merged_object_vim_pos <- FindClusters(merged_object_vim_pos, resolution = 0.05)		
merged_object_vim_pos <- RunUMAP(merged_object_vim_pos, dims = 1:35)
DimPlot(merged_object_vim_pos, reduction = "umap")
DimPlot(merged_object_vim_pos, reduction = "umap", group.by = c("orig.ident", "seurat_clusters"), split.by = "orig.ident")

# proliferative cell markers --------------------------------------------
plot <- FeaturePlot(merged_object_vim_pos, features = c("Ube2c"))
ggsave(filename = '/Users/kodeebao/Documents/Final Analysis v_2/COLD and CON cells, Vim+/35 dims, 0.05 res/proliferative_feature.png', plot, width = 8, height = 5)

plot <- VlnPlot(merged_object_vim_pos, features = c("Ube2c"))
ggsave(filename = '/Users/kodeebao/Documents/Final Analysis v_2/COLD and CON cells, Vim+/35 dims, 0.05 res/proliferative_vln.png', plot, width = 8, height = 5)


# RGCs-glia cell markers --------------------------------------------
plot <- FeaturePlot(merged_object_vim_pos, features = c("Gdpd2"))
ggsave(filename = '/Users/kodeebao/Documents/Final Analysis v_2/COLD and CON cells, Vim+/35 dims, 0.05 res/RGCs-glia_feature.png', plot, width = 8, height = 5)

plot <- VlnPlot(merged_object_vim_pos, features = c("Gdpd2"))
ggsave(filename = '/Users/kodeebao/Documents/Final Analysis v_2/COLD and CON cells, Vim+/35 dims, 0.05 res/RGCs-glia_vln.png', plot, width = 8, height = 5)


# OPCs cell markers --------------------------------------------
plot <- FeaturePlot(merged_object_vim_pos, features = c("Olig2",
                                                "Sox10",
                                                "Pdgfra"))
ggsave(filename = '/Users/kodeebao/Documents/Final Analysis v_2/COLD and CON cells, Vim+/35 dims, 0.05 res/OPCs_feature.png', plot, width = 20, height = 10)

plot <- VlnPlot(merged_object_vim_pos, features = c("Olig2",
                                            "Sox10",
                                            "Pdgfra"))
ggsave(filename = '/Users/kodeebao/Documents/Final Analysis v_2/COLD and CON cells, Vim+/35 dims, 0.05 res/OPCs_vln.png', plot, width = 20, height = 10)


# Astrocytes cell markers --------------------------------------------
plot <- FeaturePlot(merged_object_vim_pos, features = c("Aldh1l1",
                                                "Gfap",
                                                "Slc1a3"))
ggsave(filename = '/Users/kodeebao/Documents/Final Analysis v_2/COLD and CON cells, Vim+/35 dims, 0.05 res/astrocytes_feature.png', plot, width = 20, height = 10)

plot <- VlnPlot(merged_object_vim_pos, features = c("Aldh1l1",
                                            "Gfap",
                                            "Slc1a3"))
ggsave(filename = '/Users/kodeebao/Documents/Final Analysis v_2/COLD and CON cells, Vim+/35 dims, 0.05 res/astrocytes_vln.png', plot, width = 20, height = 10)

# blackshaw NPC markers cell markers --------------------------------------------
plot <- FeaturePlot(merged_object_vim_pos, features = c("Ascl1",
                                                "Neurog2"))
ggsave(filename = '/Users/kodeebao/Documents/Final Analysis v_2/COLD and CON cells, Vim+/35 dims, 0.05 res/blackshaw_NPCs_feature.png', plot, width = 15, height = 5)

plot <- VlnPlot(merged_object_vim_pos, features = c("Ascl1",
                                            "Neurog2"))
ggsave(filename = '/Users/kodeebao/Documents/Final Analysis v_2/COLD and CON cells, Vim+/35 dims, 0.05 res/blackshaw_NPCs_vln.png', plot, width = 15, height = 5)


# ependymal markers cell markers --------------------------------------------
plot <- FeaturePlot(merged_object_vim_pos, features = c("Rarres2",
                                                "Foxj1"))
ggsave(filename = '/Users/kodeebao/Documents/Final Analysis v_2/COLD and CON cells, Vim+/35 dims, 0.05 res/ependymal_feature.png', plot, width = 15, height = 5)

plot <- VlnPlot(merged_object_vim_pos, features = c("Rarres2",
                                            "Foxj1"))
ggsave(filename = '/Users/kodeebao/Documents/Final Analysis v_2/COLD and CON cells, Vim+/35 dims, 0.05 res/ependymal_vln.png', plot, width = 15, height = 5)


# tanycytes markers cell markers --------------------------------------------
plot <- FeaturePlot(merged_object_vim_pos, features = c("Lhx2",
                                                "Rax"))
ggsave(filename = '/Users/kodeebao/Documents/Final Analysis v_2/COLD and CON cells, Vim+/35 dims, 0.05 res/tanycytes_feature.png', plot, width = 15, height = 5)

plot <- VlnPlot(merged_object_vim_pos, features = c("Lhx2",
                                            "Rax"))
ggsave(filename = '/Users/kodeebao/Documents/Final Analysis v_2/COLD and CON cells, Vim+/35 dims, 0.05 res/tanycytes_vln.png', plot, width = 15, height = 5)



# testing harmonys markers on 35 dims, 0.1 res (merged and Vim+ subsetted) ------------------------

merged_object_vim_pos <- FindNeighbors(merged_object_vim_pos, dims = 1:35)	
merged_object_vim_pos <- FindClusters(merged_object_vim_pos, resolution = 0.1)		
merged_object_vim_pos <- RunUMAP(merged_object_vim_pos, dims = 1:35)
DimPlot(merged_object_vim_pos, reduction = "umap")
DimPlot(merged_object_vim_pos, reduction = "umap", group.by = c("orig.ident", "seurat_clusters"), split.by = "orig.ident")

# proliferative cell markers --------------------------------------------
plot <- FeaturePlot(merged_object_vim_pos, features = c("Ube2c"))
ggsave(filename = '/Users/kodeebao/Documents/Final Analysis v_2/COLD and CON cells, Vim+/35 dims, 0.1 res/proliferative_feature.png', plot, width = 8, height = 5)

plot <- VlnPlot(merged_object_vim_pos, features = c("Ube2c"))
ggsave(filename = '/Users/kodeebao/Documents/Final Analysis v_2/COLD and CON cells, Vim+/35 dims, 0.1 res/proliferative_vln.png', plot, width = 8, height = 5)


# RGCs-glia cell markers --------------------------------------------
plot <- FeaturePlot(merged_object_vim_pos, features = c("Gdpd2"))
ggsave(filename = '/Users/kodeebao/Documents/Final Analysis v_2/COLD and CON cells, Vim+/35 dims, 0.1 res/RGCs-glia_feature.png', plot, width = 8, height = 5)

plot <- VlnPlot(merged_object_vim_pos, features = c("Gdpd2"))
ggsave(filename = '/Users/kodeebao/Documents/Final Analysis v_2/COLD and CON cells, Vim+/35 dims, 0.1 res/RGCs-glia_vln.png', plot, width = 8, height = 5)


# OPCs cell markers --------------------------------------------
plot <- FeaturePlot(merged_object_vim_pos, features = c("Olig2",
                                                        "Sox10",
                                                        "Pdgfra"))
ggsave(filename = '/Users/kodeebao/Documents/Final Analysis v_2/COLD and CON cells, Vim+/35 dims, 0.1 res/OPCs_feature.png', plot, width = 20, height = 10)

plot <- VlnPlot(merged_object_vim_pos, features = c("Olig2",
                                                    "Sox10",
                                                    "Pdgfra"))
ggsave(filename = '/Users/kodeebao/Documents/Final Analysis v_2/COLD and CON cells, Vim+/35 dims, 0.1 res/OPCs_vln.png', plot, width = 20, height = 10)


# Astrocytes cell markers --------------------------------------------
plot <- FeaturePlot(merged_object_vim_pos, features = c("Aldh1l1",
                                                        "Gfap",
                                                        "Slc1a3"))
ggsave(filename = '/Users/kodeebao/Documents/Final Analysis v_2/COLD and CON cells, Vim+/35 dims, 0.1 res/astrocytes_feature.png', plot, width = 20, height = 10)

plot <- VlnPlot(merged_object_vim_pos, features = c("Aldh1l1",
                                                    "Gfap",
                                                    "Slc1a3"))
ggsave(filename = '/Users/kodeebao/Documents/Final Analysis v_2/COLD and CON cells, Vim+/35 dims, 0.1 res/astrocytes_vln.png', plot, width = 20, height = 10)

# blackshaw NPC markers cell markers --------------------------------------------
plot <- FeaturePlot(merged_object_vim_pos, features = c("Ascl1",
                                                        "Neurog2"))
ggsave(filename = '/Users/kodeebao/Documents/Final Analysis v_2/COLD and CON cells, Vim+/35 dims, 0.1 res/blackshaw_NPCs_feature.png', plot, width = 15, height = 5)

plot <- VlnPlot(merged_object_vim_pos, features = c("Ascl1",
                                                    "Neurog2"))
ggsave(filename = '/Users/kodeebao/Documents/Final Analysis v_2/COLD and CON cells, Vim+/35 dims, 0.1 res/blackshaw_NPCs_vln.png', plot, width = 15, height = 5)


# ependymal markers cell markers --------------------------------------------
plot <- FeaturePlot(merged_object_vim_pos, features = c("Rarres2",
                                                        "Foxj1"))
ggsave(filename = '/Users/kodeebao/Documents/Final Analysis v_2/COLD and CON cells, Vim+/35 dims, 0.1 res/ependymal_feature.png', plot, width = 15, height = 5)

plot <- VlnPlot(merged_object_vim_pos, features = c("Rarres2",
                                                    "Foxj1"))
ggsave(filename = '/Users/kodeebao/Documents/Final Analysis v_2/COLD and CON cells, Vim+/35 dims, 0.1 res/ependymal_vln.png', plot, width = 15, height = 5)


# tanycytes markers cell markers --------------------------------------------
plot <- FeaturePlot(merged_object_vim_pos, features = c("Lhx2",
                                                        "Rax"))
ggsave(filename = '/Users/kodeebao/Documents/Final Analysis v_2/COLD and CON cells, Vim+/35 dims, 0.1 res/tanycytes_feature.png', plot, width = 15, height = 5)

plot <- VlnPlot(merged_object_vim_pos, features = c("Lhx2",
                                                    "Rax"))
ggsave(filename = '/Users/kodeebao/Documents/Final Analysis v_2/COLD and CON cells, Vim+/35 dims, 0.1 res/tanycytes_vln.png', plot, width = 15, height = 5)

