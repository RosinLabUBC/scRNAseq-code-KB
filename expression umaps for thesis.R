library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)


# done on all subsetted cells (CON & COLD, M & F)

# load in object 
clean_combined_object <- readRDS('/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/clean_combined_object_CON&COLD_strat2.rds')

# Cycling
plot <- FeaturePlot(clean_combined_object, features = c("Birc5", 
                                                        "Ccnb1",
                                                        "Sgo2a",
                                                        "Cdca8",
                                                        "Tpx2",
                                                        "Ube2c"))
ggsave(filename = '/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Vim+ and subsetted cells/expression umaps of all cells for thesis/cycling/cycling_feature4.png', plot, width = 10, height = 12)


# NPCs-Glia
plot <- FeaturePlot(clean_combined_object, features = c("Cntn1",
                                                        "Olig1",
                                                        "Olig2",
                                                        "Pdgfra",
                                                        "C1ql1",
                                                        "Lhfpl3"))
ggsave(filename = '/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Vim+ and subsetted cells/expression umaps of all cells for thesis/npcs-glia/npcs-glia_feature3.png', plot, width = 10, height = 12)


# NPCs-Neuronal
plot <- FeaturePlot(clean_combined_object, features = c("Dlx1",
                                                        "Dlx2",
                                                        "Dlx5",
                                                        "Arx",
                                                        "Sp9",
                                                        "Lhx6"))
ggsave(filename = '/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Vim+ and subsetted cells/expression umaps of all cells for thesis/npcs-neuronal/npcs-neuronal_feature3.png', plot, width = 10, height = 12)


# Glial Precursors
plot <- FeaturePlot(clean_combined_object, features = c("Aldh1l1",
                                                        "Aldoc",
                                                        "Fabp7",
                                                        "Hes5",
                                                        "Mfge8",
                                                        "Slc1a3"))
ggsave(filename = '/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Vim+ and subsetted cells/expression umaps of all cells for thesis/glial precursors/glial precursors_feature4.png', plot, width = 10, height = 12)


# Neuronal Precursors (Ascl1+)
plot <- FeaturePlot(clean_combined_object, features = c("Ascl1",
                                                        "Miat",
                                                        "Isl1",
                                                        "Pak3",
                                                        "Magel2",
                                                        "Tub"))
ggsave(filename = '/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Vim+ and subsetted cells/expression umaps of all cells for thesis/ascl1+/ascl1_feature2.png', plot, width = 10, height = 12)


# Neuronal Precursors (Neurog2+) 
plot <- FeaturePlot(clean_combined_object, features = c("Lmx1a",
                                                        "Foxa1",
                                                        "Neurod1",
                                                        "Neurod2",
                                                        "Neurog2",
                                                        "Nhlh1"))
ggsave(filename = '/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Vim+ and subsetted cells/expression umaps of all cells for thesis/neurog2+/neurog2_feature2.png', plot, width = 10, height = 12)


# testing
plot <- FeaturePlot(clean_combined_object, features = c("Gad1",
                                                        "Dlx1",
                                                        "Dlx2",
                                                        "Dlx5",
                                                        "Lhx6"))
ggsave(filename = '/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Vim+ and subsetted cells/expression umaps of all cells for thesis/gad1.png', plot, width = 10, height = 12)

