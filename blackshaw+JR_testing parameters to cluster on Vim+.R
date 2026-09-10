# testing harmonys markers on 35 dims, 0.1 res (merged and Vim+ subsetted) ------------------------

merged_object_vim_pos <- FindNeighbors(merged_object_vim_pos, dims = 1:35)	
merged_object_vim_pos <- FindClusters(merged_object_vim_pos, resolution = 0.1)		
merged_object_vim_pos <- RunUMAP(merged_object_vim_pos, dims = 1:35)
DimPlot(merged_object_vim_pos, reduction = "umap")
DimPlot(merged_object_vim_pos, reduction = "umap", group.by = c("orig.ident", "seurat_clusters"), split.by = "orig.ident")

# proliferative cell markers --------------------------------------------
plot <- FeaturePlot(merged_object_vim_pos, features = c("Ube2c"))
ggsave(filename = '/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/testing parameters on Vim+/35 dims, 0.1 res/proliferative_feature.png', plot, width = 8, height = 5)

plot <- VlnPlot(merged_object_vim_pos, features = c("Ube2c"))
ggsave(filename = '/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/testing parameters on Vim+/35 dims, 0.1 res/proliferative_vln.png', plot, width = 8, height = 5)


# RGCs-glia cell markers --------------------------------------------
plot <- FeaturePlot(merged_object_vim_pos, features = c("Gdpd2"))
ggsave(filename = '/Users/kodeebao/Documents/Final Analysis v_2/COLD and CON cells, Vim+/35 dims, 0.1 res/RGCs-glia_feature.png', plot, width = 8, height = 5)

plot <- VlnPlot(merged_object_vim_pos, features = c("Gdpd2"))
ggsave(filename = '/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/testing parameters on Vim+/35 dims, 0.1 res/RGCs-glia_vln.png', plot, width = 8, height = 5)


# OPCs cell markers --------------------------------------------
plot <- FeaturePlot(merged_object_vim_pos, features = c("Olig2",
                                                        "Sox10",
                                                        "Pdgfra"))
ggsave(filename = '/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/testing parameters on Vim+/35 dims, 0.1 res/OPCs_feature.png', plot, width = 20, height = 10)

plot <- VlnPlot(merged_object_vim_pos, features = c("Olig2",
                                                    "Sox10",
                                                    "Pdgfra"))
ggsave(filename = '/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/testing parameters on Vim+/35 dims, 0.1 res/OPCs_vln.png', plot, width = 20, height = 10)


# Astrocytes cell markers --------------------------------------------
plot <- FeaturePlot(merged_object_vim_pos, features = c("Aldh1l1",
                                                        "Gfap",
                                                        "Slc1a3"))
ggsave(filename = '/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/testing parameters on Vim+/35 dims, 0.1 res/astrocytes_feature.png', plot, width = 20, height = 10)

plot <- VlnPlot(merged_object_vim_pos, features = c("Aldh1l1",
                                                    "Gfap",
                                                    "Slc1a3"))
ggsave(filename = '/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/testing parameters on Vim+/35 dims, 0.1 res/astrocytes_vln.png', plot, width = 20, height = 10)

# blackshaw NPC markers cell markers --------------------------------------------
plot <- FeaturePlot(merged_object_vim_pos, features = c("Ascl1",
                                                        "Neurog2"))
ggsave(filename = '/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/testing parameters on Vim+/35 dims, 0.1 res/blackshaw_NPCs_feature.png', plot, width = 15, height = 5)

plot <- VlnPlot(merged_object_vim_pos, features = c("Ascl1",
                                                    "Neurog2"))
ggsave(filename = '/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/testing parameters on Vim+/35 dims, 0.1 res/blackshaw_NPCs_vln.png', plot, width = 15, height = 5)


# ependymal markers cell markers --------------------------------------------
plot <- FeaturePlot(merged_object_vim_pos, features = c("Rarres2",
                                                        "Foxj1"))
ggsave(filename = '/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/testing parameters on Vim+/35 dims, 0.1 res/ependymal_feature.png', plot, width = 15, height = 5)

plot <- VlnPlot(merged_object_vim_pos, features = c("Rarres2",
                                                    "Foxj1"))
ggsave(filename = '/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/testing parameters on Vim+/35 dims, 0.1 res/ependymal_vln.png', plot, width = 15, height = 5)


# tanycytes markers cell markers --------------------------------------------
plot <- FeaturePlot(merged_object_vim_pos, features = c("Lhx2",
                                                        "Rax"))
ggsave(filename = '/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/testing parameters on Vim+/35 dims, 0.1 res/tanycytes_feature.png', plot, width = 15, height = 5)

plot <- VlnPlot(merged_object_vim_pos, features = c("Lhx2",
                                                    "Rax"))
ggsave(filename = '/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/testing parameters on Vim+/35 dims, 0.1 res/tanycytes_vln.png', plot, width = 15, height = 5)


# testing harmonys markers on 35 dims, 0.05 res (merged and Vim+ subsetted) ------------------------

merged_object_vim_pos <- FindNeighbors(merged_object_vim_pos, dims = 1:35)	
merged_object_vim_pos <- FindClusters(merged_object_vim_pos, resolution = 0.05)		
merged_object_vim_pos <- RunUMAP(merged_object_vim_pos, dims = 1:35)
DimPlot(merged_object_vim_pos, reduction = "umap")
DimPlot(merged_object_vim_pos, reduction = "umap", group.by = c("orig.ident", "seurat_clusters"), split.by = "orig.ident")

# proliferative cell markers --------------------------------------------
plot <- FeaturePlot(merged_object_vim_pos, features = c("Ube2c"))
ggsave(filename = '/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/testing parameters on Vim+/35 dims, 0.05 res/proliferative_feature.png', plot, width = 8, height = 5)

plot <- VlnPlot(merged_object_vim_pos, features = c("Ube2c"))
ggsave(filename = '/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/testing parameters on Vim+/35 dims, 0.05 res/proliferative_vln.png', plot, width = 8, height = 5)


# RGCs-glia cell markers --------------------------------------------
plot <- FeaturePlot(merged_object_vim_pos, features = c("Gdpd2"))
ggsave(filename = '/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/testing parameters on Vim+/35 dims, 0.05 res/RGCs-glia_feature.png', plot, width = 8, height = 5)

plot <- VlnPlot(merged_object_vim_pos, features = c("Gdpd2"))
ggsave(filename = '/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/testing parameters on Vim+/35 dims, 0.05 res/RGCs-glia_vln.png', plot, width = 8, height = 5)


# OPCs cell markers --------------------------------------------
plot <- FeaturePlot(merged_object_vim_pos, features = c("Olig2",
                                                        "Sox10",
                                                        "Pdgfra"))
ggsave(filename = '/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/testing parameters on Vim+/35 dims, 0.05 res/OPCs_feature.png', plot, width = 20, height = 10)

plot <- VlnPlot(merged_object_vim_pos, features = c("Olig2",
                                                    "Sox10",
                                                    "Pdgfra"))
ggsave(filename = '/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/testing parameters on Vim+/35 dims, 0.05 res/OPCs_vln.png', plot, width = 20, height = 10)


# Astrocytes cell markers --------------------------------------------
plot <- FeaturePlot(merged_object_vim_pos, features = c("Aldh1l1",
                                                        "Gfap",
                                                        "Slc1a3"))
ggsave(filename = '/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/testing parameters on Vim+/35 dims, 0.05 res/astrocytes_feature.png', plot, width = 20, height = 10)

plot <- VlnPlot(merged_object_vim_pos, features = c("Aldh1l1",
                                                    "Gfap",
                                                    "Slc1a3"))
ggsave(filename = '/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/testing parameters on Vim+/35 dims, 0.05 res/astrocytes_vln.png', plot, width = 20, height = 10)

# blackshaw NPC markers cell markers --------------------------------------------
plot <- FeaturePlot(merged_object_vim_pos, features = c("Ascl1",
                                                        "Neurog2"))
ggsave(filename = '/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/testing parameters on Vim+/35 dims, 0.05 res/blackshaw_NPCs_feature.png', plot, width = 15, height = 5)

plot <- VlnPlot(merged_object_vim_pos, features = c("Ascl1",
                                                    "Neurog2"))
ggsave(filename = '/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/testing parameters on Vim+/35 dims, 0.05 res/blackshaw_NPCs_vln.png', plot, width = 15, height = 5)


# ependymal markers cell markers --------------------------------------------
plot <- FeaturePlot(merged_object_vim_pos, features = c("Rarres2",
                                                        "Foxj1"))
ggsave(filename = '/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/testing parameters on Vim+/35 dims, 0.05 res/ependymal_feature.png', plot, width = 15, height = 5)

plot <- VlnPlot(merged_object_vim_pos, features = c("Rarres2",
                                                    "Foxj1"))
ggsave(filename = '/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/testing parameters on Vim+/35 dims, 0.05 res/ependymal_vln.png', plot, width = 15, height = 5)


# tanycytes markers cell markers --------------------------------------------
plot <- FeaturePlot(merged_object_vim_pos, features = c("Lhx2",
                                                        "Rax"))
ggsave(filename = '/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/testing parameters on Vim+/35 dims, 0.05 res/tanycytes_feature.png', plot, width = 15, height = 5)

plot <- VlnPlot(merged_object_vim_pos, features = c("Lhx2",
                                                    "Rax"))
ggsave(filename = '/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/testing parameters on Vim+/35 dims, 0.05 res/tanycytes_vln.png', plot, width = 15, height = 5)


# testing harmonys markers on 33 dims, 0.05 res (merged and Vim+ subsetted) ------------------------

merged_object_vim_pos <- FindNeighbors(merged_object_vim_pos, dims = 1:33)	
merged_object_vim_pos <- FindClusters(merged_object_vim_pos, resolution = 0.05)		
merged_object_vim_pos <- RunUMAP(merged_object_vim_pos, dims = 1:33)
DimPlot(merged_object_vim_pos, reduction = "umap")
DimPlot(merged_object_vim_pos, reduction = "umap", group.by = c("orig.ident", "seurat_clusters"), split.by = "orig.ident")

# proliferative cell markers --------------------------------------------
plot <- FeaturePlot(merged_object_vim_pos, features = c("Ube2c"))
ggsave(filename = '/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/testing parameters on Vim+/33 dims, 0.05 res/proliferative_feature.png', plot, width = 8, height = 5)

plot <- VlnPlot(merged_object_vim_pos, features = c("Ube2c"))
ggsave(filename = '/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/testing parameters on Vim+/33 dims, 0.05 res/proliferative_vln.png', plot, width = 8, height = 5)


# RGCs-glia cell markers --------------------------------------------
plot <- FeaturePlot(merged_object_vim_pos, features = c("Gdpd2"))
ggsave(filename = '/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/testing parameters on Vim+/33 dims, 0.05 res/RGCs-glia_feature.png', plot, width = 8, height = 5)

plot <- VlnPlot(merged_object_vim_pos, features = c("Gdpd2"))
ggsave(filename = '/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/testing parameters on Vim+/33 dims, 0.05 res/RGCs-glia_vln.png', plot, width = 8, height = 5)


# OPCs cell markers --------------------------------------------
plot <- FeaturePlot(merged_object_vim_pos, features = c("Olig2",
                                                        "Sox10",
                                                        "Pdgfra"))
ggsave(filename = '/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/testing parameters on Vim+/33 dims, 0.05 res/OPCs_feature.png', plot, width = 20, height = 10)

plot <- VlnPlot(merged_object_vim_pos, features = c("Olig2",
                                                    "Sox10",
                                                    "Pdgfra"))
ggsave(filename = '/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/testing parameters on Vim+/33 dims, 0.05 res/OPCs_vln.png', plot, width = 20, height = 10)


# Astrocytes cell markers --------------------------------------------
plot <- FeaturePlot(merged_object_vim_pos, features = c("Aldh1l1",
                                                        "Gfap",
                                                        "Slc1a3"))
ggsave(filename = '/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/testing parameters on Vim+/33 dims, 0.05 res/astrocytes_feature.png', plot, width = 20, height = 10)

plot <- VlnPlot(merged_object_vim_pos, features = c("Aldh1l1",
                                                    "Gfap",
                                                    "Slc1a3"))
ggsave(filename = '/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/testing parameters on Vim+/33 dims, 0.05 res/astrocytes_vln.png', plot, width = 20, height = 10)

# blackshaw NPC markers cell markers --------------------------------------------
plot <- FeaturePlot(merged_object_vim_pos, features = c("Ascl1",
                                                        "Neurog2"))
ggsave(filename = '/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/testing parameters on Vim+/33 dims, 0.05 res/blackshaw_NPCs_feature.png', plot, width = 15, height = 5)

plot <- VlnPlot(merged_object_vim_pos, features = c("Ascl1",
                                                    "Neurog2"))
ggsave(filename = '/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/testing parameters on Vim+/33 dims, 0.05 res/blackshaw_NPCs_vln.png', plot, width = 15, height = 5)


# ependymal markers cell markers --------------------------------------------
plot <- FeaturePlot(merged_object_vim_pos, features = c("Rarres2",
                                                        "Foxj1"))
ggsave(filename = '/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/testing parameters on Vim+/33 dims, 0.05 res/ependymal_feature.png', plot, width = 15, height = 5)

plot <- VlnPlot(merged_object_vim_pos, features = c("Rarres2",
                                                    "Foxj1"))
ggsave(filename = '/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/testing parameters on Vim+/33 dims, 0.05 res/ependymal_vln.png', plot, width = 15, height = 5)


# tanycytes markers cell markers --------------------------------------------
plot <- FeaturePlot(merged_object_vim_pos, features = c("Lhx2",
                                                        "Rax"))
ggsave(filename = '/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/testing parameters on Vim+/33 dims, 0.05 res/tanycytes_feature.png', plot, width = 15, height = 5)

plot <- VlnPlot(merged_object_vim_pos, features = c("Lhx2",
                                                    "Rax"))
ggsave(filename = '/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/testing parameters on Vim+/33 dims, 0.05 res/tanycytes_vln.png', plot, width = 15, height = 5)

