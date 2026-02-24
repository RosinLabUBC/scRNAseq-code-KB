# using bubble plot to graph treatment groups against eachother rather than individually 

library(CellChat)
library(patchwork)
library(ggplot2)
library(dplyr)
library(Seurat)
library(future)

# Load preprocessed CellChat objects
female_con_NSCclusters_mg_cellchat <- readRDS("/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/CellChat Objects/NSC clusters vs microglia/female_con_NSCclusters_mg_cellchat.rds")
female_cold_NSCclusters_mg_cellchat <- readRDS("/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/CellChat Objects/NSC clusters vs microglia/female_cold_NSCclusters_mg_cellchat.rds")
male_con_NSCclusters_mg_cellchat <- readRDS("/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/CellChat Objects/NSC clusters vs microglia/male_con_NSCclusters_mg_cellchat.rds")
male_cold_NSCclusters_mg_cellchat <- readRDS("/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/CellChat Objects/NSC clusters vs microglia/male_cold_NSCclusters_mg_cellchat.rds")


# con vs cold based on sex ------------------------------------------------
# Male: CON vs COLD
male.list <- list(
  CON  = male_con_NSCclusters_mg_cellchat,
  COLD = male_cold_NSCclusters_mg_cellchat)

male.merge <- mergeCellChat(
  male.list,
  add.names = names(male.list))

# Female: CON vs COLD
female.list <- list(
  CON  = female_con_NSCclusters_mg_cellchat,
  COLD = female_cold_NSCclusters_mg_cellchat)

female.merge <- mergeCellChat(
  female.list,
  add.names = names(female.list))

#cycling to cycling
netVisual_bubble(
  male.merge,
  sources.use  = "Cycling NSPCs",
  targets.use  = "Cycling NSPCs",
  comparison   = c(1, 2),   # CON vs COLD
  remove.isolate = TRUE,
  angle.x = 45,
  title.name = "Male: CON vs COLD (Cycling NSPCs → Cycling NSPCs)")

netVisual_bubble(
  female.merge,
  sources.use  = "Cycling NSPCs",
  targets.use  = "Cycling NSPCs",
  comparison   = c(1, 2),
  remove.isolate = TRUE,
  angle.x = 45,
  title.name = "Female: CON vs COLD (Cycling NSPCs → Cycling NSPCs)")

#microglia to cyclling
netVisual_bubble(
  male.merge,
  sources.use = "Microglia",
  targets.use = "Cycling NSPCs",
  comparison = c(1,2),
  remove.isolate = TRUE,
  angle.x = 45,
  title.name = "Male: CON vs COLD (Microglia → Cycling NSPCs)")

netVisual_bubble(
  female.merge,
  sources.use = "Microglia",
  targets.use = "Cycling NSPCs",
  comparison = c(1,2),
  remove.isolate = TRUE,
  angle.x = 45,
  title.name = "Female: CON vs COLD (Microglia → Cycling NSPCs)")

#cycling to microglia
netVisual_bubble(
  male.merge,
  sources.use = "Cycling NSPCs",
  targets.use = "Microglia",
  comparison = c(1,2),
  remove.isolate = TRUE,
  angle.x = 45,
  title.name = "Male: CON vs COLD (Cycling NSPCs → Microglia)")

netVisual_bubble(
  female.merge,
  sources.use = "Cycling NSPCs",
  targets.use = "Microglia",
  comparison = c(1,2),
  remove.isolate = TRUE,
  angle.x = 45,
  title.name = "Female: CON vs COLD (Cycling NSPCs → Microglia)")


# all 4 sexes and conditions together for thesis -------------------------------------

all.list <- list(
  Male_CON   = male_con_NSCclusters_mg_cellchat,
  Male_COLD  = male_cold_NSCclusters_mg_cellchat,
  Female_CON = female_con_NSCclusters_mg_cellchat,
  Female_COLD= female_cold_NSCclusters_mg_cellchat)

all.merge <- mergeCellChat(
  all.list,
  add.names = names(all.list))

# cycling to cycling
netVisual_bubble(
  all.merge,
  sources.use  = "Cycling NSPCs",
  targets.use  = "Cycling NSPCs",
  comparison   = c(1, 2, 3, 4),  # Male_CON, Male_COLD, Female_CON, Female_COLD
  remove.isolate = TRUE,
  angle.x = 45,
  title.name = "Cycling NSPCs → Cycling NSPCs (Male & Female, CON vs COLD)")

# microglia to cycling
netVisual_bubble(
  all.merge,
  sources.use  = "Microglia",
  targets.use  = "Cycling NSPCs",
  comparison   = c(1, 2, 3, 4),
  remove.isolate = TRUE,
  angle.x = 45,
  title.name = "Microglia → Cycling NSPCs (Male & Female, CON vs COLD)")

# microglia to cycling
netVisual_bubble(
  all.merge,
  sources.use  = "Cycling NSPCs",
  targets.use  = "Microglia",
  comparison   = c(1, 2, 3, 4),
  remove.isolate = TRUE,
  angle.x = 45,
  title.name = "Cycling NSPCs → Microglia (Male & Female, CON vs COLD)")


