library(dplyr)
library(tidyverse)
library(IOBR)
library(immunedeconv)
library(RColorBrewer)
library(ggrepel)
library(clusterProfiler)
library(msigdbr)
library(AnnotationDbi)

options(ggrepel.max.overlaps = Inf)
set.seed(1)

# set working directory
wd <- "~/B_Memory_master/Master_Diploma_Private/Output/Bulk//"
setwd(wd)

# Reading in the feature count file 
feature_counts <- read.table("Data/featurecounts_results.txt", 
                             header = TRUE, row.names = 1) %>%
  select(6:ncol(.))

# Renaming the column names
colnames(feature_counts) <- sub("_alignedAligned.sortedByCoord.out.bam", "", 
                            sub("X.home.eshchoka.B_Memory_master.Master_Diploma_Private.data.alignment2.bam.", "",
                            sub("IgG1", "IgG_1", colnames(feature_counts))))

# Creating metadata dataframe
metadata <- data.frame(file = colnames(feature_counts)) %>%
  separate(file, into = c("cancer_type", "patient", "tissue_type", "isotype"))
rownames(metadata) <- colnames(feature_counts)


# Selecting Rcp_tum
feature_counts_Rcp_tum <- feature_counts %>% 
  select(contains("Rcp")) %>% 
  select(contains("tum"))

metadata_Rcp_tum <- metadata %>% 
  filter(cancer_type == "Rcp") %>% 
  select(-cancer_type) %>% 
  filter(tissue_type == "tum") %>% 
  select(-tissue_type)

# Deconvolution -----------------------------------------------------------

# TPM_normalization
feature_counts_Rcp_tum_tpm <- count2tpm(feature_counts_Rcp_tum,
                                        idType = "Ensembl",
                                        org = "hsa")

# Deconvolution with xCell
feature_counts_Rcp_tum_xCell <- deconvolute(feature_counts_Rcp_tum_tpm, "xcell") %>% 
  as.data.frame() %>%
  tibble::column_to_rownames('cell_type')

# Creating a heatmap
pdf("Graphs/Rcp_xCell.pdf")
pheatmap(as.data.frame(feature_counts_Rcp_tum_xCell),
         cluster_rows = T,
         show_rownames = T,
         border_color = NA,
         fontsize = 10,
         fontsize_row = 8,
         height = 20)
dev.off()

# Correlation between MS4A1 and CD19
df_Rcp <- data.frame("MS4A1" = as.numeric(feature_counts_Rcp_tum_tpm["MS4A1",]),
                 "CD19" = as.numeric(feature_counts_Rcp_tum_tpm["CD19",]),
                 row.names = colnames(feature_counts_Rcp_tum_tpm))

# Plotting correlation
pdf("Graphs/Rcp_CD19_vs_CD20.pdf")
ggplot(df_Rcp, aes(MS4A1, CD19)) + 
  geom_point() + 
  geom_text_repel(aes(label=rownames(df_Rcp))) +
  geom_smooth(method=lm, se=FALSE, size=0.5) +
  theme_bw()
dev.off()

# Filtering samples
samples_to_remove <- c("Rcp_15_tum_IgA_2_S23",
                       "Rcp_15_tum_IgG_2_S24",
                       "Rcp_8_tum_IgG_2_S65", # 34.9% assigned in featureCounts
                       "Rcp19_tum_IgG_S5", 
                       "Rcp19_tum_IgA_S6", 
                       "Rcp20_tum_IgG_S7", 
                       "Rcp20_tum_IgA_S8",
                       "Rcp16_tum_IgG1_S50",
                       "Rcp_16_tum_IgA_S51",
                       "Rcp_17_tum_IgG_S41",
                       "Rcp_12_tum_IgG_S37",
                       "Rcp_22_tum_IgG_1_S11",
                       "Rcp_22_tum_IgG_2_S12",
                       "Rcp_22_tum_IgA_1_S13",
                       "Rcp_22_tum_IgA_2_S14",
                       "Rcp_18_tum_IgA_1_star_S3",
                       "Rcp_18_tum_IgG_1_star_S1",
                       "Rcp_18_tum_IgA_2_S4")

feature_counts_Rcp_tum_final <- select(feature_counts_Rcp_tum, -any_of(samples_to_remove))
metadata_Rcp_tum_final <- metadata_Rcp_tum[rownames(metadata_Rcp_tum) %in% colnames(feature_counts_Rcp_tum_final),]


# DE analysis -------------------------------------------------------------

# Begin Differential Expression (DE) analysis
dds_Rcp_tum <- DESeqDataSetFromMatrix(countData = feature_counts_Rcp_tum_final, 
                                      colData = metadata_Rcp_tum_final, 
                                      design = ~ patient + isotype)
dds_Rcp_tum <- estimateSizeFactors(dds_Rcp_tum)

# Filter out genes where there are less than 3 samples with normalized counts greater than or equal to 5
idx_Rcp_tum <- rowSums(counts(dds_Rcp_tum, normalized=TRUE) >= 5 ) >= 3 
dds_Rcp_tum <- dds_Rcp_tum[idx_Rcp_tum,]

dds_Rcp_tum <- DESeq(dds_Rcp_tum)
rld_Rcp_tum <- vst(dds_Rcp_tum, blind=F)

# Plot PCA plot
pdf("Graphs/Rcp_PCA.pdf")
plotPCA(rld_Rcp_tum, intgroup="isotype") + 
  geom_text_repel(aes(label=name)) + 
  ggtitle("Rcp_tum") + 
  theme_bw()
dev.off()

summary(results(dds_Rcp_tum, alpha=0.05)) #

contrast <- c("isotype", "IgA", "IgG")
res_unshrunken_Rcp_tum <- results(dds_Rcp_tum, contrast=contrast)

# Plot MA before shrinkage
pdf("Graphs/Rcp_LFC_before_shrinkage.pdf")
plotMA(res_unshrunken_Rcp_tum, main="No shrinkage of LFCs", ylim=c(-3,3), cex=1)
dev.off()

res_Rcp_tum <- lfcShrink(dds_Rcp_tum, coef="isotype_IgG_vs_IgA", type="apeglm")

# Plot MA after shrinkage
pdf("Graphs/Rcp_LFC_after_apeglm_shrinkage.pdf")
plotMA(res_Rcp_tum, main="Shrinkage of LFCs(apeglm)", ylim=c(-3,3), cex=1)
dev.off()

# Convert results to a data frame and map gene names
res_table_Rcp_tum <- res_Rcp_tum %>% 
  as.data.frame() %>% 
  mutate(refs = mapIds(org.Hs.eg.db, keys = rownames(.), column = c('SYMBOL'), keytype = 'ENSEMBL')) %>% 
  drop_na()

# Prepare data for the volcano plot

# Filter data from uninformative genes
res_table_Rcp_tum <- res_table_Rcp_tum[- grep("IGH|IGK|IGL|LOC|JCHAIN|LINC", res_table_Rcp_tum$refs),]

# Add a column of NAs
res_table_Rcp_tum$diffexpressed <- "Unchanged"
# Set as "UP" if log2Foldchange > 1 and padj < 0.05
res_table_Rcp_tum$diffexpressed[res_table_Rcp_tum$log2FoldChange > 1 & res_table_Rcp_tum$padj < 0.05] <- "Up-regulated"
# Set as "DOWN" if log2Foldchange < -1 and padj < 0.05
res_table_Rcp_tum$diffexpressed[res_table_Rcp_tum$log2FoldChange < -1 & res_table_Rcp_tum$padj < 0.05] <- "Down-regulated"



# Create volcano plot and save as PDF

# Select top up- and down-regulated genes
top_n = 20
top_genes_Rcp_tum <- bind_rows(res_table_Rcp_tum %>% 
                         dplyr::filter(diffexpressed == 'Up-regulated') %>%  
                         dplyr::top_n(-top_n, wt = padj) %>% 
                         head(top_n),
                       res_table_Rcp_tum %>% 
                         dplyr::filter(diffexpressed == 'Down-regulated') %>%  
                         dplyr::top_n(-top_n, wt = padj) %>% 
                         head(top_n))

# Create the volcano plot
p_Rcp <- ggplot(res_table_Rcp_tum) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = diffexpressed), size = 7/5) +
  scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
  ggtitle("A vs G(p_adj<0.05 & |LFC|>1)") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") + 
  geom_vline(xintercept=c(-1, 1), col="black", linetype="dotted") +
  geom_hline(yintercept=-log10(0.05), col="black", linetype="dotted") +
  guides(colour = guide_legend(override.aes = list(size=1.5))) +
  geom_label_repel(data = top_genes_Rcp_tum,
                   mapping = aes(log2FoldChange, -log(padj,10), label = refs),
                   size = 3.5) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme_bw()

# Save the plot to a Pres_table_Lcp_tum file
ggsave("Graphs/Rcp_volcano_plot.pdf", plot = p_Rcp, height = 6, width = 10)


# Filter the significant results for only upregulated and downregulated genes
significant_res_Rcp_tum <- dplyr::filter(res_table_Rcp_tum, padj < 0.05 & abs(log2FoldChange) > 1)

# Create a new object with the log-transformed normalized count data for the significant genes
rld_counts_sig_Rcp_tum <- assay(rld_Rcp_tum)[rownames(significant_res_Rcp_tum), ] %>%
  as.data.frame() %>% 
  mutate(refs = mapIds(org.Hs.eg.db, keys = rownames(.), column = c('SYMBOL'), keytype = 'ENSEMBL'))%>% 
  drop_na() %>% 
  remove_rownames() %>% 
  tibble::column_to_rownames("refs")

# Create a heatmap of the significant genes
pdf("Graphs/Rcp_heatmap.pdf")
pheatmap(
  rld_counts_sig_Rcp_tum,
  cluster_rows = T,
  show_rownames = T,
  annotation_col = metadata_Rcp_tum_final,
  border_color = "grey20",
  fontsize = 10,
  scale = "row",
  fontsize_row = 8
)
dev.off()

# Order the significant results by adjusted p-value
significant_res_Rcp_tum <- significant_res_Rcp_tum[order(significant_res_Rcp_tum$log2FoldChange), ] %>%
  remove_rownames() %>% 
  tibble::column_to_rownames("refs")

# Export the results to a TSV file
write.table(
  significant_res_Rcp_tum,
  file = "Tables/Rcp_DE_genes.tsv",
  sep = "\t",
  quote = FALSE
)


# GSEA --------------------------------------------------------------------


hs_hallmark_sets <- msigdbr(
  species = "Homo sapiens", 
  category = "H"
)

if(any(duplicated(res_table_Rcp_tum$refs))) {
  
  # Arrange by absolute log fold change in descending order and then by p-value in ascending order
  res_table_Rcp_tum <- res_table_Rcp_tum %>%
    dplyr::arrange(-abs(log2FoldChange), pvalue)
  
  # Retain only the first occurrence of each ref (which will be the one with the highest absolute LFC and lowest p-value)
  res_table_Rcp_tum <- res_table_Rcp_tum %>%
    dplyr::distinct(refs, .keep_all = TRUE)
}

# Let's create a named vector ranked based on the log2 fold change values
lfc_vector <- res_table_Rcp_tum$log2FoldChange
names(lfc_vector) <- res_table_Rcp_tum$refs

# We need to sort the log2 fold change values in descending order here
lfc_vector <- sort(lfc_vector, decreasing = TRUE)

# Use the combined score for GSEA
gsea_results <- GSEA(
  geneList = lfc_vector, # Use the combined score vector
  minGSSize = 25,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  eps = 0,
  seed = TRUE,
  pAdjustMethod = "BH",
  TERM2GENE = dplyr::select(hs_hallmark_sets, gs_name, gene_symbol)
)

# Write to .tsv files
write.table(gsea_results@result, 
            file = "Tables/Rcp_gsea_results.tsv", 
            sep = "\t", 
            row.names = T, 
            quote = FALSE)

# Run GSEA using the GO database
gsea_go_BP_results <- gseGO(
  geneList = lfc_vector,  
  OrgDb = org.Hs.eg.db,
  keyType = 'SYMBOL',
  ont = "BP", 
  minGSSize = 10, 
  maxGSSize = 500, 
  pvalueCutoff = 0.05, 
  verbose = FALSE
)
gsea_go_MF_results <- gseGO(
  geneList = lfc_vector,  
  OrgDb = org.Hs.eg.db,
  keyType = 'SYMBOL',
  ont = "MF", 
  minGSSize = 10, 
  maxGSSize = 500, 
  pvalueCutoff = 0.05, 
  verbose = FALSE
)
gsea_go_CC_results <- gseGO(
  geneList = lfc_vector,  
  OrgDb = org.Hs.eg.db,
  keyType = 'SYMBOL',
  ont = "CC", 
  minGSSize = 10, 
  maxGSSize = 500, 
  pvalueCutoff = 0.05, 
  verbose = FALSE
)

# Write to .tsv files
write.table(gsea_go_BP_results@result, 
            file = "Tables/Rcp_gsea_go_BP_results.tsv", 
            sep = "\t", 
            row.names = T, 
            quote = FALSE)
write.table(gsea_go_MF_results@result, 
            file = "Tables/Rcp_gsea_go_MF_results.tsv", 
            sep = "\t", 
            row.names = T, 
            quote = FALSE)
write.table(gsea_go_CC_results@result, 
            file = "Tables/Rcp_gsea_go_CC_results.tsv", 
            sep = "\t", 
            row.names = T, 
            quote = FALSE)

# Visualize GSEA output
for (variable in gsea_results@result$ID) {
  tmp <- enrichplot::gseaplot(
    gsea_results,
    geneSetID = variable,
    title = gsub("^HALLMARK_", "", variable),
    color.line = "#0d76ff"
  )
  
  # Modify the PDF file name with the specific index
  pdf_file <- paste("Graphs/Rcp_GSEA_", gsub("^HALLMARK_", "", variable), ".pdf", sep = "")
  
  # Save the plot to the PDF file using ggsave()
  ggsave(pdf_file, plot = tmp, height = 14, width = 10)
  
}

# Visualize GO output
output_list <- c(gsea_go_BP_results, gsea_go_MF_results, gsea_go_CC_results)

variable <- gsea_go_BP_results
for (variable in output_list) {
  for (counter in 1:nrow(variable@result)) {
    tmp <- enrichplot::gseaplot(
      variable,
      geneSetID = variable@result[counter,1],
      title = variable@result[counter,2],
      color.line = "#0d76ff"
    )
    
    # Modify the PDF file name with the specific index
    pdf_file <- paste("Graphs/Rcp_GSEA_GO_", variable@setType ,"_" , gsub(" ", "_", variable@result[counter,2]), ".pdf", sep = "")
    
    # Save the plot to the PDF file using ggsave()
    ggsave(pdf_file, plot = tmp, height = 14, width = 10)
  }
}