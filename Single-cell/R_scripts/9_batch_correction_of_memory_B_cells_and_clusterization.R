library(Seurat)
library(tidyverse)
library(patchwork)
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

# Helper function for feature plotting with custom colors
featurePlotCustom <- function(data, feature) {
  FeaturePlot(data, features = feature, reduction = "umap", pt.size = 0.05, raster=FALSE) + 
    scale_colour_gradientn(colours = c("grey80", "#3366CC", "#0000FF", "#0000CC", "#000099"))
}

# Function to save plot to png
saveToPng <- function(plot, path_suffix) {
  png(paste0("Graphs_png/9_", path_suffix, ".png"), 
      width = 10, 
      height = 7, 
      res = 450, 
      units = "in")
  print(plot)
  dev.off()
}

# Function to save plot to png at lower resolution
saveToPng_small <- function(plot, path_suffix) {
  png(paste0("Graphs_png/9_", path_suffix, "_small.png"), 
      width = 6, 
      height = 4, 
      res = 450, 
      units = "in")
  print(plot)
  dev.off()
}

load("8_Corrected_B_cells.RData")

int_2 <- subset(immune.combined, idents = c(5, 6, 7, 11), invert = TRUE)

DefaultAssay(int_2) <- "RNA"

int_2.list <- SplitObject(int_2, split.by = "patient")
rm(int_2)

# select var feat
data(cc.genes)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
var.genes_pattern <- "TR[ABGD][CV]|RP[SL]|IG[HKL][VDJ]|IGHM|IGHD|IGHE|IGHA[1-2]|IGHG[1-4]|IGKC|IGLC_mem[1-7]|AC233755.1|IGLL|JCHAIN|MT-|MALAT1|MTRNR|XIST|IGHGP|HSP|DNAJ|JUN|DUSP|IER|FOS|BTG|ATF3|NR4A1"
var.genes <- grep(pattern = var.genes_pattern, x = rownames(int_2.list), value = TRUE)
excl_genes <- c(s.genes, g2m.genes, var.genes)

# normalize and identify variable features for each dataset independently
int_2.list <- lapply(X = int_2.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 1000)
  VariableFeatures(x) <- VariableFeatures(x)[-which(VariableFeatures(x) %in% excl_genes)]
  x
})

# select features that are repeatedly variable across datasets for integration
int_2.features <- SelectIntegrationFeatures(object.list = int_2.list)
int_2.anchors <- FindIntegrationAnchors(object.list = int_2.list, anchor.features = int_2.features, normalization.method = "LogNormalize")
rm(int_2.list)

# this command creates an 'integrated' data assay
int_2 <- IntegrateData(anchorset = int_2.anchors, normalization.method = "LogNormalize", k.weight = 15)

# Switch to integrated assay.
DefaultAssay(int_2) <- "integrated"

int_2 <- ScaleData(int_2, verbose = FALSE)
int_2 <- RunPCA(int_2, verbose = FALSE)
int_2 <- RunUMAP(int_2, reduction = "pca", dims = 1:15)
int_2 <- FindNeighbors(int_2, reduction = "pca", dims = 1:15)
int_2 <- FindClusters(int_2, resolution = 1.2)


saveToPng(DimPlot(int_2, reduction = "umap", group.by = "patient"), "DimPlot_patient")
saveToPng_small(DimPlot(int_2, reduction = "umap", group.by = "patient"), "DimPlot_patient")

saveToPng(DimPlot(int_2, reduction = "umap", label = TRUE, repel = TRUE), "DimPlot_clusters")
saveToPng_small(DimPlot(int_2, reduction = "umap", label = TRUE, repel = TRUE), "DimPlot_clusters")

DefaultAssay(int_2) <- "RNA"

Names <- c("FCRL4", "PDCD1", "CR2", "TNFSF11", "FAS", "ITGAX", "CR2", "TBX21")

# Saving FeaturePlots to png
lapply(Names, function(feat) {
  plot <- featurePlotCustom(int_2, feat)
  saveToPng(plot, paste0("FeaturePlot_", feat))
  saveToPng_small(plot, paste0("FeaturePlot_", feat))
})

# Saving VlnPlots to png
lapply(Names, function(feat) {
  plot <- VlnPlot(int_2, feat) +
    stat_summary(fun = mean, geom='point', size = 5, colour = "black", shape = 95)
  saveToPng(plot, paste0("VlnPlot_", feat))
  saveToPng_small(plot, paste0("VlnPlot_", feat))
})

# search for DE genes among all clusters
all_markers <- FindAllMarkers(int_2) %>% filter(p_val_adj < 0.001)
write.table(all_markers[, c(ncol(all_markers), 1:(ncol(all_markers)-1))], "Tables/9_all_DE_genes.tsv", row.names = F, col.names = T, sep = "\t")

# rm all accept int_2
rm(list = ls()[-match("int_2", ls())])

# Save workspace
save.image("9_Corrected_memory_B_cells.RData")
