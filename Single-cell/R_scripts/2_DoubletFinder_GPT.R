# Libraries
library(Seurat)
library(DoubletFinder)
library(tidyverse)

options(ggrepel.max.overlaps = Inf)
set.seed(1)

# set working directory
wd <- "~/B_Memory_master/Master_Diploma_Private/Output/Single-cell/"
setwd(wd)

# Load Data
load("1_Labeled.RData")

# Split the object into samples
LC.split <- SplitObject(LC, split.by = "orig.ident")
rm(LC)

# Loop through samples to find doublets
for(i in 1:length(LC.split)) {
  LC.sample <- LC.split[[i]] %>%
    NormalizeData() %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA(nfeatures.print = 10)
  
  # Find significant PCs
  stdv <- LC.sample[["pca"]]@stdev
  sum.stdv <- sum(stdv)
  percent.stdv <- (stdv / sum.stdv) * 100
  cumulative <- cumsum(percent.stdv)
  co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
  co2 <- sort(which((percent.stdv[1:length(percent.stdv) - 1] -
                       percent.stdv[2:length(percent.stdv)]) > 0.1),
              decreasing = T)[1] + 1
  min.pc <- min(co1, co2)
  
  # UMAP, Neighbors, and Clusters
  LC.sample <- LC.sample %>%
    RunUMAP(dims = 1:min.pc) %>%
    FindNeighbors(dims = 1:min.pc) %>%
    FindClusters(resolution = 0.1)
  
  # Doublet Finding
  sweep.list <- paramSweep_v3(LC.sample, PCs = 1:min.pc, num.cores = detectCores() - 1)
  sweep.stats <- summarizeSweep(sweep.list)
  bcmvn <- find.pK(sweep.stats)

  bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric),]
  optimal.pk <- bcmvn.max$pK
  optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk] # Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
  
  annotations <- LC.sample@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations) ## Homotypic doublet proportion estimate
  nExp.poi <- round(optimal.pk * nrow(LC.sample@meta.data)) ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))
  
  LC.sample <- doubletFinder_v3(seu = LC.sample,
                                PCs = 1:min.pc,
                                pK = optimal.pk,
                                nExp = nExp.poi.adj)
  
  # Finalize Sample
  colnames(LC.sample@meta.data)[11] <- "doublet_finder"
  LC.singlets <- LC.sample
  LC.singlets@meta.data <- LC.singlets@meta.data[c(1:7, 11)] # to remove empty metadata columns
  LC.split[[i]] <- LC.singlets
  rm(list = ls()[-match("LC.split", ls())])
}

# Merge Samples
LC <- LC.split[[1]]
for (i in 2:length(LC.split)) {
  LC <- merge(x = LC, y = LC.split[[i]])
}

# rm all accept LC
rm(list = ls()[-match("LC", ls())])

# save workspace
save.image("2_Labeled_doublets.RData")
