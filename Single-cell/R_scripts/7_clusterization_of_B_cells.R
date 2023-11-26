# Libraries
library(Seurat)
library(ggplot2)
library(dplyr)
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

# Subset data for memory cells
LC_B <- subset(LC, idents = c(1, 25, 52, 65))

# Clear unnecessary variables
rm(list = setdiff(ls(), "LC_B"))

# Utility function to save plots to png
save_plot <- function(plot_obj, filename) {
  png(paste0("Graphs_png/7_", filename, ".png"), 
      width = 10, 
      height = 7, 
      res = 450, 
      units = "in")
  print(plot_obj)
  dev.off()
}

# Analysis
LC_B <- FindVariableFeatures(LC_B, selection.method = "vst", nfeatures = 1000)

# Show and save most variable features
a <- LabelPoints(plot = VariableFeaturePlot(LC_B), points = head(VariableFeatures(LC_B), 30), repel = FALSE)
save_plot(a, "Dotplot_variable_features")

# Filter genes using dplyr
data(cc.genes)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

var.genes_pattern <- "TR[ABGD][CV]|RP[SL]|IG[HKL][VDJ]|IGHM|IGHD|IGHE|IGHA[1-2]|IGHG[1-4]|IGKC|IGLC_B[1-7]|AC233755.1|IGLL|JCHAIN|MT-|MALAT1|MTRNR|XIST|IGHGP|HSP|DNAJ|JUN|DUSP|IER|FOS|BTG|ATF3|NR4A1"
var.genes <- grep(pattern = var.genes_pattern, x = rownames(LC_B), value = TRUE)

VariableFeatures(LC_B) %>% 
  setdiff(var.genes) %>% 
  setdiff(s.genes) %>% 
  setdiff(g2m.genes) -> VariableFeatures(LC_B)

# Show and save most variable features after filtration
b <- LabelPoints(plot = VariableFeaturePlot(LC_B), points = head(VariableFeatures(LC_B), 30), repel = FALSE)
save_plot(b, "Dotplot_variable_features_after_filtration")

LC_B <- ScaleData(LC_B, features = rownames(LC_B))
LC_B <- RunPCA(LC_B, features = VariableFeatures(object = LC_B))
LC_B <- FindNeighbors(LC_B, dims = 1:20)
LC_B <- FindClusters(LC_B, resolution = 1)
LC_B <- RunUMAP(LC_B, dims = 1:20)

# Save plots
getPalette_clust = colorRampPalette(brewer.pal(9, "Set1"))(length(unique(LC_B$seurat_clusters)))
c <- DimPlot(LC_B, reduction = "umap", label = TRUE, cols = getPalette_clust, raster=FALSE) + NoLegend()
save_plot(c, "Dim_plot")
getPalette_clust_pat = colorRampPalette(brewer.pal(9, "Set1"))(length(unique(LC_B$patient)))
d <- DimPlot(LC_B, reduction = "umap", group.by = "patient", cols = getPalette_clust_pat, raster=FALSE)
save_plot(d, "Dim_plot_patients")

# rm all accept LC_B
rm(list = ls()[-match("LC_B", ls())])

# Save workspace
save.image("7_Clustered_B_cells.RData")
