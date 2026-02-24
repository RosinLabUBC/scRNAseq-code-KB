# combined object after the overlapping and unassigned cells were removed

library(SeuratExtend)
library(Seurat)
library(dplyr)
library(patchwork)
library(clustree)
library(ggplot2)

# for thesis figure ---------------------------------------------

# load in object
clean_combined_object <- readRDS('/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/clean_combined_object_CON&COLD_strat2.rds')

# Add combined group column
clean_combined_object$group <- paste(clean_combined_object$sex, clean_combined_object$treatment, sep = "_")

# Set factor levels for consistent group ordering
clean_combined_object$group <- factor(
  clean_combined_object$group,
  levels = c("Male_CON", "Male_COLD", "Female_CON", "Female_COLD"))

# Define custom cluster colors
cluster_colors <- c(
  "Neuronal Precursors (Ascl1+)" = "#F8766D",
  "Glial Precursors" = "#B79F00",
  "Cycling NSPCs" = "#00BA38",
  "NPCs-Glia" = "#00BFC4",
  "NPCs-Neuronal" = "#649DFF",
  "Neuronal Precursors (Neurog2+)" = "#F564E3")

clean_combined_object$cluster_names <- factor(
  clean_combined_object$cluster_names,
  levels = c( "Cycling NSPCs",
              "NPCs-Glia",
              "NPCs-Neuronal",
              "Glial Precursors",
              "Neuronal Precursors (Ascl1+)",
              "Neuronal Precursors (Neurog2+)"))

# Plot cluster distribution barplot
plot <- ClusterDistrBar(
  origin = clean_combined_object$group,
  cluster = clean_combined_object$cluster_names,
  rev = FALSE,
  normalize = TRUE,
  percent = TRUE,
  plot = TRUE,
  flip = TRUE,
  reverse_order = TRUE)

# Apply cluster colors
plot + scale_fill_manual(values = cluster_colors)

