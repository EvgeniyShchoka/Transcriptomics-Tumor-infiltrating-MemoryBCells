# Libraries
library(Seurat)
library(RColorBrewer)

options(ggrepel.max.overlaps = Inf)
set.seed(1)

# set working directory
wd <- "~/B_Memory_master/Master_Diploma_Private/Output/Single-cell/"
setwd(wd)

# Create directory for output if it doesn't exist
if (!dir.exists("Graphs_png")) {
  dir.create("Graphs_png")
}

load("5_Clustered_all_cells.RData") 

# Helper function to get color palette
getPalette <- function(data_col, n) {
  len <- length(unique(data_col))
  return(colorRampPalette(brewer.pal(min(n, 9), "Set1"))(len))
}

# Helper function for feature plotting with custom colors
featurePlotCustom <- function(data, feature) {
  FeaturePlot(data, features = feature, reduction = "umap", pt.size = 0.05, raster=FALSE) + 
    scale_colour_gradientn(colours = c("grey80", "#3366CC", "#0000FF", "#0000CC", "#000099"))
}

# Function to save plot to png
saveToPng <- function(plot, path_suffix) {
  png(paste0("Graphs_png/6_", path_suffix, ".png"), 
      width = 10, 
      height = 7, 
      res = 450, 
      units = "in")
  print(plot)
  dev.off()
}

# Function to save plot to png at lower resolution
saveToPng_small <- function(plot, path_suffix) {
  png(paste0("Graphs_png/6_", path_suffix, "_small.png"), 
      width = 6, 
      height = 4, 
      res = 450, 
      units = "in")
  print(plot)
  dev.off()
}

# Saving DimPlots to png
plots <- list(
  DimPlot(LC, reduction = "umap", label = TRUE, cols = getPalette(LC$seurat_clusters, 9), raster=FALSE) + NoLegend(),
  DimPlot(LC, reduction = "umap", label = TRUE, group.by = "patient", cols = getPalette(LC$patient, 9), raster=FALSE),
  DimPlot(LC, reduction = "umap", label = TRUE, group.by = "lib", cols = getPalette(LC$lib, 3), raster=FALSE),
  DimPlot(LC, reduction = "umap", group.by = "doublet_finder", cols = getPalette(LC$doublet_finder, 3), label = FALSE, raster=FALSE),
  DimPlot(LC, reduction = "umap", label = TRUE, group.by = "cancer_type", cols = getPalette(LC$cancer_type, 3), raster=FALSE)
)
paths <- c("clusters", "patient", "lib", "doublet_finder", "cancer_type")

for(i in seq_along(paths)) {
  saveToPng(plots[[i]], paste0("DimPlot_", paths[i]))
  saveToPng_small(plots[[i]], paste0("DimPlot_", paths[i]))
}

# Saving FeaturePlots to png
featureNames <- c("MS4A1", "CD3E", "CD38", "CD19", "FOS")
lapply(featureNames, function(feat) {
  plot <- featurePlotCustom(LC, feat)
  saveToPng(plot, paste0("FeaturePlot_", feat))
  saveToPng_small(plot, paste0("FeaturePlot_", feat))
})
