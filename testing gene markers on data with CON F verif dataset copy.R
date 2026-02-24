
# Ascl1+ ------------------------------------------------------------------

# proliferative cell markers
plot <- FeaturePlot(merged_object_vim_pos, features = c("Ascl1"))
ggsave(filename = '/Users/kodeebao/Documents/Strategy 2 Final Analysis/CON F Verification 2/Vim+/testing gene markers_all vim+/35 dims, 0.1 res/ascl1_feature.png', plot, width = 8, height = 5)

plot <- VlnPlot(merged_object_vim_pos, features = c("Ascl1"))
ggsave(filename = '/Users/kodeebao/Documents/Strategy 2 Final Analysis/CON F Verification 2/Vim+/testing gene markers_all vim+/35 dims, 0.1 res/ascl1_vln.png', plot, width = 8, height = 5)

# Seths Markers -----------------------------------------------------------

# proliferative cell markers
plot <- FeaturePlot(merged_object_vim_pos, features = c("Ube2c"))
ggsave(filename = '/Users/kodeebao/Documents/Strategy 2 Final Analysis/CON F Verification 2/Vim+/testing gene markers_all vim+/35 dims, 0.1 res/Seths Markers_proliferative_feature.png', plot, width = 8, height = 5)

plot <- VlnPlot(merged_object_vim_pos, features = c("Ube2c"))
ggsave(filename = '/Users/kodeebao/Documents/Strategy 2 Final Analysis/CON F Verification 2/Vim+/testing gene markers_all vim+/35 dims, 0.1 res/Seths Markers_proliferative_vln.png', plot, width = 8, height = 5)


# Harmony Old Markers -----------------------------------------------------

# proliferative cell markers

plot <- FeaturePlot(merged_object_vim_pos, features = c("Pif1",
                                                        "Depdc1a",
                                                        "Troap",
                                                        "Cep55",
                                                        "Aspm",
                                                        "Hmmr",
                                                        "Esco2",
                                                        "Cdc25c",
                                                        "Plk1",
                                                        "Psrc1",
                                                        "Hist1h3c",
                                                        "Ccnb1",
                                                        "Kif2c",
                                                        "Ccna2",
                                                        "Cdc20",
                                                        "Sgo1",
                                                        "Nusap1",
                                                        "Mis18bp1",
                                                        "Ube2c",
                                                        "Tpx2",
                                                        "Spc25",
                                                        "Kif18a",
                                                        "Prc1",
                                                        "Hist1h3b",
                                                        "Sgo2a",
                                                        "Kifc1",
                                                        "Cenpf",
                                                        "Bub1b",
                                                        "Dlgap5",
                                                        "Kif23"))
ggsave(filename = '/Users/kodeebao/Documents/Strategy 2 Final Analysis/CON F Verification 2/Vim+/testing gene markers_all vim+/35 dims, 0.1 res/Harmony Old_proliferative_feature.png', plot, width = 40, height = 30)

plot <- VlnPlot(merged_object_vim_pos, features = c("Pif1",
                                                    "Depdc1a",
                                                    "Troap",
                                                    "Cep55",
                                                    "Aspm",
                                                    "Hmmr",
                                                    "Esco2",
                                                    "Cdc25c",
                                                    "Plk1",
                                                    "Psrc1",
                                                    "Hist1h3c",
                                                    "Ccnb1",
                                                    "Kif2c",
                                                    "Ccna2",
                                                    "Cdc20",
                                                    "Sgo1",
                                                    "Nusap1",
                                                    "Mis18bp1",
                                                    "Ube2c",
                                                    "Tpx2",
                                                    "Spc25",
                                                    "Kif18a",
                                                    "Prc1",
                                                    "Hist1h3b",
                                                    "Sgo2a",
                                                    "Kifc1",
                                                    "Cenpf",
                                                    "Bub1b",
                                                    "Dlgap5",
                                                    "Kif23"))
ggsave(filename = '/Users/kodeebao/Documents/Strategy 2 Final Analysis/CON F Verification 2/Vim+/testing gene markers_all vim+/35 dims, 0.1 res/Harmony Old_proliferative_vln.png', plot, width = 40, height = 30)

# Harmony New Markers (no p-value cutoff) ---------------------------------
# new markers (don't need p-value < 0.05)

# Cycling RGCs 1 cell markers

plot <- FeaturePlot(merged_object_vim_pos, features = c("Cdc20",
                                                        "Ccnb2",
                                                        "Hmmr",
                                                        "Cenpa",
                                                        "Cdkn3",
                                                        "Psrc1",
                                                        "Cdc25c",
                                                        "Aspm",
                                                        "Kif2c",
                                                        "Ccnb1",
                                                        "Ube2c",
                                                        "Cep55",
                                                        "Pimreg",
                                                        "Kif20a",
                                                        "Cenpf",
                                                        "Plk1",
                                                        "Knstrn",
                                                        "Tpx2",
                                                        "Tuba1c",
                                                        "Cenpe",
                                                        "Dlgap5",
                                                        "Birc5",
                                                        "Depdc1a",
                                                        "Cdca8",
                                                        "Bub1b",
                                                        "Cdca3",
                                                        "Ect2",
                                                        "Ccna2",
                                                        "Sgo2a",
                                                        "Aurka"))
ggsave(filename = '/Users/kodeebao/Documents/Strategy 2 Final Analysis/CON F Verification 2/Vim+/testing gene markers_all vim+/35 dims, 0.1 res/Harmony New (no p-value cutoff)_Cycling_RGCs_1_feature.png', plot, width = 40, height = 30)

plot <- VlnPlot(merged_object_vim_pos, features = c("Cdc20",
                                                    "Ccnb2",
                                                    "Hmmr",
                                                    "Cenpa",
                                                    "Cdkn3",
                                                    "Psrc1",
                                                    "Cdc25c",
                                                    "Aspm",
                                                    "Kif2c",
                                                    "Ccnb1",
                                                    "Ube2c",
                                                    "Cep55",
                                                    "Pimreg",
                                                    "Kif20a",
                                                    "Cenpf",
                                                    "Plk1",
                                                    "Knstrn",
                                                    "Tpx2",
                                                    "Tuba1c",
                                                    "Cenpe",
                                                    "Dlgap5",
                                                    "Birc5",
                                                    "Depdc1a",
                                                    "Cdca8",
                                                    "Bub1b",
                                                    "Cdca3",
                                                    "Ect2",
                                                    "Ccna2",
                                                    "Sgo2a",
                                                    "Aurka"))
ggsave(filename = '/Users/kodeebao/Documents/Strategy 2 Final Analysis/CON F Verification 2/Vim+/testing gene markers_all vim+/35 dims, 0.1 res/Harmony New (no p-value cutoff)_Cycling_RGCs_1_vln.png', plot, width = 40, height = 30)


# Cycling RGCs 2 cell markers

plot <- FeaturePlot(merged_object_vim_pos, features = c("Hist1h1b",
                                                        "Esco2",
                                                        "Hist1h3c",
                                                        "Hist1h2ap",
                                                        "Tk1",
                                                        "Hist1h2ad",
                                                        "Hist1h1a",
                                                        "Neil3",
                                                        "Cdca5",
                                                        "Lrr1",
                                                        "Exo1",
                                                        "Pclaf",
                                                        "Clspn",
                                                        "Dscc1",
                                                        "Aurkb",
                                                        "Tcf19",
                                                        "Hist1h2ae",
                                                        "Rad51ap1",
                                                        "Ccne2",
                                                        "Fbxo5",
                                                        "Hist1h1d",
                                                        "Hist1h2bm",
                                                        "Rrm2",
                                                        "Gmnn",
                                                        "Asf1b",
                                                        "Ttk",
                                                        "Ndc80",
                                                        "Shcbp1",
                                                        "Ska1",
                                                        "Spc24"))
ggsave(filename = '/Users/kodeebao/Documents/Strategy 2 Final Analysis/CON F Verification 2/Vim+/testing gene markers_all vim+/35 dims, 0.1 res/Harmony New (no p-value cutoff)_Cycling_RGCs_2_feature.png', plot, width = 40, height = 30)

plot <- VlnPlot(merged_object_vim_pos, features = c("Hist1h1b",
                                                    "Esco2",
                                                    "Hist1h3c",
                                                    "Hist1h2ap",
                                                    "Tk1",
                                                    "Hist1h2ad",
                                                    "Hist1h1a",
                                                    "Neil3",
                                                    "Cdca5",
                                                    "Lrr1",
                                                    "Exo1",
                                                    "Pclaf",
                                                    "Clspn",
                                                    "Dscc1",
                                                    "Aurkb",
                                                    "Tcf19",
                                                    "Hist1h2ae",
                                                    "Rad51ap1",
                                                    "Ccne2",
                                                    "Fbxo5",
                                                    "Hist1h1d",
                                                    "Hist1h2bm",
                                                    "Rrm2",
                                                    "Gmnn",
                                                    "Asf1b",
                                                    "Ttk",
                                                    "Ndc80",
                                                    "Shcbp1",
                                                    "Ska1",
                                                    "Spc24"))
ggsave(filename = '/Users/kodeebao/Documents/Strategy 2 Final Analysis/CON F Verification 2/Vim+/testing gene markers_all vim+/35 dims, 0.1 res/Harmony New (no p-value cutoff)_Cycling_RGCs_2_vln.png', plot, width = 40, height = 30)





# Seths Markers -----------------------------------------------------------

# proliferative cell markers
plot <- FeaturePlot(merged_object_vim_pos_subset, features = c("Ube2c"))
ggsave(filename = '/Users/kodeebao/Documents/Strategy 2 Final Analysis/CON F Verification 2/Vim+/testing gene markers_vim+ subset/30 dims, 0.18 res/Seths Markers_proliferative_feature.png', plot, width = 8, height = 5)

plot <- VlnPlot(merged_object_vim_pos_subset, features = c("Ube2c"))
ggsave(filename = '/Users/kodeebao/Documents/Strategy 2 Final Analysis/CON F Verification 2/Vim+/testing gene markers_vim+ subset/30 dims, 0.18 res/Seths Markers_proliferative_vln.png', plot, width = 8, height = 5)


# Harmony Old Markers -----------------------------------------------------

# proliferative cell markers

plot <- FeaturePlot(merged_object_vim_pos_subset, features = c("Pif1",
                                                               "Depdc1a",
                                                               "Troap",
                                                               "Cep55",
                                                               "Aspm",
                                                               "Hmmr",
                                                               "Esco2",
                                                               "Cdc25c",
                                                               "Plk1",
                                                               "Psrc1",
                                                               "Hist1h3c",
                                                               "Ccnb1",
                                                               "Kif2c",
                                                               "Ccna2",
                                                               "Cdc20",
                                                               "Sgo1",
                                                               "Nusap1",
                                                               "Mis18bp1",
                                                               "Ube2c",
                                                               "Tpx2",
                                                               "Spc25",
                                                               "Kif18a",
                                                               "Prc1",
                                                               "Hist1h3b",
                                                               "Sgo2a",
                                                               "Kifc1",
                                                               "Cenpf",
                                                               "Bub1b",
                                                               "Dlgap5",
                                                               "Kif23"))
ggsave(filename = '/Users/kodeebao/Documents/Strategy 2 Final Analysis/CON F Verification 2/Vim+/testing gene markers_vim+ subset/30 dims, 0.18 res/Harmony Old_proliferative_feature.png', plot, width = 40, height = 30)

plot <- VlnPlot(merged_object_vim_pos_subset, features = c("Pif1",
                                                           "Depdc1a",
                                                           "Troap",
                                                           "Cep55",
                                                           "Aspm",
                                                           "Hmmr",
                                                           "Esco2",
                                                           "Cdc25c",
                                                           "Plk1",
                                                           "Psrc1",
                                                           "Hist1h3c",
                                                           "Ccnb1",
                                                           "Kif2c",
                                                           "Ccna2",
                                                           "Cdc20",
                                                           "Sgo1",
                                                           "Nusap1",
                                                           "Mis18bp1",
                                                           "Ube2c",
                                                           "Tpx2",
                                                           "Spc25",
                                                           "Kif18a",
                                                           "Prc1",
                                                           "Hist1h3b",
                                                           "Sgo2a",
                                                           "Kifc1",
                                                           "Cenpf",
                                                           "Bub1b",
                                                           "Dlgap5",
                                                           "Kif23"))
ggsave(filename = '/Users/kodeebao/Documents/Strategy 2 Final Analysis/CON F Verification 2/Vim+/testing gene markers_vim+ subset/30 dims, 0.18 res/Harmony Old_proliferative_vln.png', plot, width = 40, height = 30)

# Harmony New Markers (no p-value cutoff) ---------------------------------
# new markers (don't need p-value < 0.05)

# Cycling RGCs 1 cell markers

plot <- FeaturePlot(merged_object_vim_pos_subset, features = c("Cdc20",
                                                               "Ccnb2",
                                                               "Hmmr",
                                                               "Cenpa",
                                                               "Cdkn3",
                                                               "Psrc1",
                                                               "Cdc25c",
                                                               "Aspm",
                                                               "Kif2c",
                                                               "Ccnb1",
                                                               "Ube2c",
                                                               "Cep55",
                                                               "Pimreg",
                                                               "Kif20a",
                                                               "Cenpf",
                                                               "Plk1",
                                                               "Knstrn",
                                                               "Tpx2",
                                                               "Tuba1c",
                                                               "Cenpe",
                                                               "Dlgap5",
                                                               "Birc5",
                                                               "Depdc1a",
                                                               "Cdca8",
                                                               "Bub1b",
                                                               "Cdca3",
                                                               "Ect2",
                                                               "Ccna2",
                                                               "Sgo2a",
                                                               "Aurka"))
ggsave(filename = '/Users/kodeebao/Documents/Strategy 2 Final Analysis/CON F Verification 2/Vim+/testing gene markers_vim+ subset/30 dims, 0.18 res/Harmony New (no p-value cutoff)_Cycling_RGCs_1_feature.png', plot, width = 40, height = 30)

plot <- VlnPlot(merged_object_vim_pos_subset, features = c("Cdc20",
                                                           "Ccnb2",
                                                           "Hmmr",
                                                           "Cenpa",
                                                           "Cdkn3",
                                                           "Psrc1",
                                                           "Cdc25c",
                                                           "Aspm",
                                                           "Kif2c",
                                                           "Ccnb1",
                                                           "Ube2c",
                                                           "Cep55",
                                                           "Pimreg",
                                                           "Kif20a",
                                                           "Cenpf",
                                                           "Plk1",
                                                           "Knstrn",
                                                           "Tpx2",
                                                           "Tuba1c",
                                                           "Cenpe",
                                                           "Dlgap5",
                                                           "Birc5",
                                                           "Depdc1a",
                                                           "Cdca8",
                                                           "Bub1b",
                                                           "Cdca3",
                                                           "Ect2",
                                                           "Ccna2",
                                                           "Sgo2a",
                                                           "Aurka"))
ggsave(filename = '/Users/kodeebao/Documents/Strategy 2 Final Analysis/CON F Verification 2/Vim+/testing gene markers_vim+ subset/30 dims, 0.18 res/Harmony New (no p-value cutoff)_Cycling_RGCs_1_vln.png', plot, width = 40, height = 30)


# Cycling RGCs 2 cell markers

plot <- FeaturePlot(merged_object_vim_pos_subset, features = c("Hist1h1b",
                                                               "Esco2",
                                                               "Hist1h3c",
                                                               "Hist1h2ap",
                                                               "Tk1",
                                                               "Hist1h2ad",
                                                               "Hist1h1a",
                                                               "Neil3",
                                                               "Cdca5",
                                                               "Lrr1",
                                                               "Exo1",
                                                               "Pclaf",
                                                               "Clspn",
                                                               "Dscc1",
                                                               "Aurkb",
                                                               "Tcf19",
                                                               "Hist1h2ae",
                                                               "Rad51ap1",
                                                               "Ccne2",
                                                               "Fbxo5",
                                                               "Hist1h1d",
                                                               "Hist1h2bm",
                                                               "Rrm2",
                                                               "Gmnn",
                                                               "Asf1b",
                                                               "Ttk",
                                                               "Ndc80",
                                                               "Shcbp1",
                                                               "Ska1",
                                                               "Spc24"))
ggsave(filename = '/Users/kodeebao/Documents/Strategy 2 Final Analysis/CON F Verification 2/Vim+/testing gene markers_vim+ subset/30 dims, 0.18 res/Harmony New (no p-value cutoff)_Cycling_RGCs_2_feature.png', plot, width = 40, height = 30)

plot <- VlnPlot(merged_object_vim_pos_subset, features = c("Hist1h1b",
                                                           "Esco2",
                                                           "Hist1h3c",
                                                           "Hist1h2ap",
                                                           "Tk1",
                                                           "Hist1h2ad",
                                                           "Hist1h1a",
                                                           "Neil3",
                                                           "Cdca5",
                                                           "Lrr1",
                                                           "Exo1",
                                                           "Pclaf",
                                                           "Clspn",
                                                           "Dscc1",
                                                           "Aurkb",
                                                           "Tcf19",
                                                           "Hist1h2ae",
                                                           "Rad51ap1",
                                                           "Ccne2",
                                                           "Fbxo5",
                                                           "Hist1h1d",
                                                           "Hist1h2bm",
                                                           "Rrm2",
                                                           "Gmnn",
                                                           "Asf1b",
                                                           "Ttk",
                                                           "Ndc80",
                                                           "Shcbp1",
                                                           "Ska1",
                                                           "Spc24"))
ggsave(filename = '/Users/kodeebao/Documents/Strategy 2 Final Analysis/CON F Verification 2/Vim+/testing gene markers_vim+ subset/30 dims, 0.18 res/Harmony New (no p-value cutoff)_Cycling_RGCs_2_vln.png', plot, width = 40, height = 30)



# Birc5 -------------------------------------------------------------------

# proliferative cell markers
plot <- FeaturePlot(merged_object_vim_pos_subset, features = c("Birc5"))
ggsave(filename = '/Users/kodeebao/Documents/Strategy 2 Final Analysis/CON F Verification 2/Vim+/testing gene markers_vim+ subset/30 dims, 0.18 res/Birc5_vim_pos_subset_feature.png', plot, width = 8, height = 5)

plot <- VlnPlot(merged_object_vim_pos_subset, features = c("Birc5"))
ggsave(filename = '/Users/kodeebao/Documents/Strategy 2 Final Analysis/CON F Verification 2/Vim+/testing gene markers_vim+ subset/30 dims, 0.18 res/Birc5_vim_pos_subset_vln.png', plot, width = 8, height = 5)

