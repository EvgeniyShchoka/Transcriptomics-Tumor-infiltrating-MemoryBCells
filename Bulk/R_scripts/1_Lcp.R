# Loading libraries
library(tidyverse)
library(pheatmap)
library(ggrepel)
library(DESeq2)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ggplot2)
library(ggpubr)
library(purrr)
library(clusterProfiler)
library(msigdbr)
library(IOBR)
library(immunedeconv)

options(ggrepel.max.overlaps = Inf)
set.seed(1)

# set working directory
wd <- "~/B_Memory_master/Master_Diploma_Private/Output/Bulk/"
setwd(wd)

# Create directory for output if it doesn't exist
if (!dir.exists("Graphs_png")) {
  dir.create("Graphs_png")
}

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

# Selecting Lcp_tum
feature_counts_Lcp_tum <- feature_counts %>% 
  select(contains("Lcp")) %>% 
  select(contains("tum"))

metadata_Lcp_tum <- metadata %>% 
  filter(cancer_type == "Lcp") %>% 
  select(-cancer_type) %>% 
  filter(tissue_type == "tum") %>% 
  select(-tissue_type)


# Deconvolution -----------------------------------------------------------


# TPM normalization
feature_counts_Lcp_tum_tpm <- count2tpm(feature_counts_Lcp_tum,
                                        idType = "Ensembl",
                                        org = "hsa")

# Deconvolution with xCell
feature_counts_Lcp_tum_xCell <- deconvolute(feature_counts_Lcp_tum_tpm, "xcell") %>% 
  as.data.frame() %>%
  tibble::column_to_rownames('cell_type')



# Creating a heatmap
png(filename = "Graphs_png/Lcp_heatmap_xCell.png", 
    width = 8, 
    height = 8, 
    res = 450, 
    units = "in")

Lcp_pheatmap_xCell <- pheatmap(feature_counts_Lcp_tum_xCell,
                               show_rownames = TRUE,
                               cluster_rows = TRUE,
                               show_row_dend = FALSE,
                               border_color = NA,
                               fontsize = 9,
                               fontsize_row = 9,
                               angle_col = "315",
                               color = colorRampPalette(c("white", "red2"))(50))

# change heading of a color palette
Lcp_pheatmap_xCell@matrix_color_mapping@name <- "cell-type score"
Lcp_pheatmap_xCell

dev.off()

# Correlation between MS4A1 and CD19
df_Lcp <- data.frame("MS4A1" = as.numeric(feature_counts_Lcp_tum_tpm["MS4A1",]),
                     "CD19" = as.numeric(feature_counts_Lcp_tum_tpm["CD19",]),
                     row.names = colnames(feature_counts_Lcp_tum_tpm))

# Plotting correlation
png(filename = "Graphs_png/Lcp_scatterplot_CD19_vs_CD20.png", 
    width = 7, 
    height = 7, 
    res = 450, 
    units = "in")

ggplot(df_Lcp, aes(MS4A1, CD19)) + 
  geom_point() + 
  geom_text_repel(aes(label=rownames(df_Lcp))) +
  geom_smooth(method=lm, se=FALSE, size=0.5) +
  theme_bw()

dev.off()

# Filtering samples
samples_to_remove <- c("Lcp_12_tum_IgA_S19",
                       "Lcp_15_tum_IgG_1_S22",
                       "Lcp_15_tum_IgA_2_S23",
                       "Lcp_15_tum_IgG_2_S24",
                       "Lcp_8_tum_IgG_2_S65")

feature_counts_Lcp_tum_final <- select(feature_counts_Lcp_tum, -any_of(samples_to_remove))
metadata_Lcp_tum_final <- metadata_Lcp_tum[rownames(metadata_Lcp_tum) %in% colnames(feature_counts_Lcp_tum_final),]


# DE analysis -------------------------------------------------------------


# Begin Differential Expression (DE) analysis
dds_Lcp_tum <- DESeqDataSetFromMatrix(countData = feature_counts_Lcp_tum_final, 
                                      colData = metadata_Lcp_tum_final, 
                                      design = ~ patient + isotype)
dds_Lcp_tum <- estimateSizeFactors(dds_Lcp_tum)

# Filter out genes where there are less than 3 samples with normalized counts greater than or equal to 5
idx_Lcp_tum <- rowSums(counts(dds_Lcp_tum, normalized=TRUE) >= 5 ) >= 3 
dds_Lcp_tum <- dds_Lcp_tum[idx_Lcp_tum,]

dds_Lcp_tum <- DESeq(dds_Lcp_tum)
rld_Lcp_tum <- vst(dds_Lcp_tum, blind=F)

# Plot PCA plot
png(filename = "Graphs_png/Lcp_PCA.png", 
    width = 7, 
    height = 7, 
    res = 450, 
    units = "in")

plotPCA(rld_Lcp_tum, intgroup="isotype") + 
  geom_text_repel(aes(label=name)) + 
  ggtitle("Lcp_tum") + 
  theme_bw()

dev.off()

summary(results(dds_Lcp_tum, alpha=0.01)) # 237

contrast <- c("isotype", "IgA", "IgG")
res_unshrunken_Lcp_tum <- results(dds_Lcp_tum, contrast=contrast)

# Plot MA before shrinkage
png(filename = "Graphs_png/Lcp_LFC_before_shrinkage.png", 
    width = 6, 
    height = 6, 
    res = 450, 
    units = "in")

plotMA(res_unshrunken_Lcp_tum, main="No shrinkage of LFCs", ylim=c(-3,3), cex=1)

dev.off()

res_Lcp_tum <- lfcShrink(dds_Lcp_tum, coef="isotype_IgG_vs_IgA", type="apeglm")

# Plot MA after shrinkage
png(filename = "Graphs_png/Lcp_LFC_after_apeglm_shrinkage.png", 
    width = 6, 
    height = 6, 
    res = 450, 
    units = "in")

plotMA(res_Lcp_tum, main="Shrinkage of LFCs(apeglm)", ylim=c(-3,3), cex=1)

dev.off()

# Invert LFC values (for further visualization)
res_Lcp_tum$log2FoldChange <- -res_Lcp_tum$log2FoldChange

# Convert results to a data frame and map gene names
res_table_Lcp_tum <- res_Lcp_tum %>% 
  as.data.frame() %>% 
  mutate(refs = mapIds(org.Hs.eg.db, keys = rownames(.), column = c('SYMBOL'), keytype = 'ENSEMBL')) %>% 
  drop_na()

# Filter data from uninformative genes
res_table_Lcp_tum <- res_table_Lcp_tum[- grep("IGH|IGK|IGL|LOC|JCHAIN|LINC", res_table_Lcp_tum$refs),]

# Add a column of NAs
res_table_Lcp_tum$diffexpressed <- "Unchanged"
# Set as "UP" if log2Foldchange > 1 and padj < 0.01
res_table_Lcp_tum$diffexpressed[res_table_Lcp_tum$log2FoldChange > 1 & res_table_Lcp_tum$padj < 0.01] <- "Up-regulated in IgA"
# Set as "DOWN" if log2Foldchange < -1 and padj < 0.01
res_table_Lcp_tum$diffexpressed[res_table_Lcp_tum$log2FoldChange < -1 & res_table_Lcp_tum$padj < 0.01] <- "Down-regulated in IgA"

# Select top up- and down-regulated genes
top_n = 20
top_genes_Lcp_tum <- bind_rows(res_table_Lcp_tum %>% 
                                 dplyr::filter(diffexpressed == 'Up-regulated in IgA') %>%  
                                 dplyr::top_n(-top_n, wt = padj) %>% 
                                 head(top_n),
                               res_table_Lcp_tum %>% 
                                 dplyr::filter(diffexpressed == 'Down-regulated in IgA') %>%  
                                 dplyr::top_n(-top_n, wt = padj) %>% 
                                 head(top_n))

# Create the volcano plot and save as png
p_Lcp <- ggplot(res_table_Lcp_tum, aes(x = log2FoldChange, y = -log10(padj), colour = diffexpressed)) +
  geom_point(size = 7/5) +
  scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
  ggtitle("IgG vs IgA (p_adj<0.01 & |LFC|>1)") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") + 
  geom_vline(xintercept=c(-1, 1), col="black", linetype="dotted") +
  geom_hline(yintercept=-log10(0.01), col="black", linetype="dotted") +
  guides(colour = guide_legend(override.aes = list(size=1.5))) +
  geom_label_repel(data = top_genes_Lcp_tum,
                   mapping = aes(log2FoldChange, -log(padj,10), label = refs),
                   size = 3.5) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme_bw()

# Save the plot
ggsave("Graphs_png/Lcp_volcano_plot_DE_genes.png", plot = p_Lcp, height = 8, width = 8)

# Filter the significant results for only upregulated and downregulated genes
significant_res_Lcp_tum <- dplyr::filter(res_table_Lcp_tum, padj < 0.01 & abs(log2FoldChange) > 1)

# Create a new object with the log-transformed normalized count data for the significant genes
rld_counts_sig_Lcp_tum <- assay(rld_Lcp_tum)[rownames(significant_res_Lcp_tum), ] %>%
  as.data.frame() %>% 
  mutate(refs = mapIds(org.Hs.eg.db, keys = rownames(.), column = c('SYMBOL'), keytype = 'ENSEMBL'))%>% 
  drop_na() %>% 
  remove_rownames() %>% 
  tibble::column_to_rownames("refs")

# Subset isotype for metadata visualization 
metadata_Lcp_tum_final_isotype <- metadata_Lcp_tum_final[, 2, drop = FALSE]

# Plot a heatmap of significant genes
png(filename = "Graphs_png/Lcp_heatmap_DE_genes.png", 
    width = 10, 
    height = 10, 
    res = 450, 
    units = "in")

Lcp_heatmap_DE_genes <- pheatmap(rld_counts_sig_Lcp_tum,
                                 color = colorRampPalette(c("navy", "white", "red"))(50),
                                 border_color = NA,
                                 cluster_rows = F,
                                 show_rownames = T,
                                 show_colnames = T,
                                 annotation_col = metadata_Lcp_tum_final_isotype,
                                 annotation_colors = list(isotype = c(IgA="forestgreen", IgG="yellow2")),
                                 fontsize = 10,
                                 fontsize_row = 10, 
                                 scale = "row",
                                 angle_col = "315")


# delete heading of a color palette
Lcp_heatmap_DE_genes@matrix_color_mapping@name <- " "
Lcp_heatmap_DE_genes

dev.off()

# Order the significant results by adjusted p-value
significant_res_Lcp_tum <- significant_res_Lcp_tum[order(significant_res_Lcp_tum$log2FoldChange), ] %>%
  remove_rownames() %>% 
  tibble::column_to_rownames("refs")

# Export the results to a TSV file
write.table(
  significant_res_Lcp_tum,
  file = "Tables/Lcp_DE_genes.tsv",
  sep = "\t",
  quote = FALSE
)


# GSEA --------------------------------------------------------------------


hs_hallmark_sets <- msigdbr(
  species = "Homo sapiens", 
  category = "H"
)

if(any(duplicated(res_table_Lcp_tum$refs))) {
  
  # Arrange by absolute log fold change in descending order and then by p-value in ascending order
  res_table_Lcp_tum <- res_table_Lcp_tum %>%
    dplyr::arrange(-abs(log2FoldChange), pvalue)
  
  # Retain only the first occurrence of each ref (which will be the one with the highest absolute LFC and lowest p-value)
  res_table_Lcp_tum <- res_table_Lcp_tum %>%
    dplyr::distinct(refs, .keep_all = TRUE)
}

# Let's create a named vector ranked based on the log2 fold change values
lfc_vector <- res_table_Lcp_tum$log2FoldChange
names(lfc_vector) <- res_table_Lcp_tum$refs

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

# 

# Write to .tsv files
write.table(gsea_results@result, 
            file = "Tables/Lcp_gsea_results.tsv", 
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
            file = "Tables/Lcp_gsea_go_BP_results.tsv", 
            sep = "\t", 
            row.names = T, 
            quote = FALSE)
write.table(gsea_go_MF_results@result, 
            file = "Tables/Lcp_gsea_go_MF_results.tsv", 
            sep = "\t", 
            row.names = T, 
            quote = FALSE)
write.table(gsea_go_CC_results@result, 
            file = "Tables/Lcp_gsea_go_CC_results.tsv", 
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
  png_file <- paste("Graphs_png/Lcp_GSEA_", gsub("^HALLMARK_", "", variable), ".png", sep = "")
  
  # Save the plot to the PDF file using ggsave()
  ggsave(png_file, plot = tmp, height = 14, width = 10)
  
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
    png_file <- paste("Graphs_png/Lcp_GSEA_GO_", variable@setType ,"_" , gsub(" ", "_", variable@result[counter,2]), ".png", sep = "")
    
    # Save the plot to the png file using ggsave()
    ggsave(png_file, plot = tmp, height = 14, width = 10)
  }
}


# Combine all the GSEA result data frames into one, with an additional 'Method' column
combined_results <- bind_rows(
  mutate(gsea_results@result, Method = "GSEA"),
  mutate(gsea_go_BP_results@result, Method = "GSEA GO BP"),
  mutate(gsea_go_MF_results@result, Method = "GSEA GO MF"),
  mutate(gsea_go_CC_results@result, Method = "GSEA GO CC")
)

# Create a list of uninformative pathways
row_names_df_to_remove<-c("HALLMARK_XENOBIOTIC_METABOLISM",
                          "HALLMARK_UV_RESPONSE_UP",
                          "HALLMARK_MYOGENESIS",
                          "HALLMARK_ESTROGEN_RESPONSE_LATE",
                          "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", 
                          "HALLMARK_COAGULATION", 
                          "HALLMARK_CHOLESTEROL_HOMEOSTASIS", 
                          "HALLMARK_BILE_ACID_METABOLISM", 
                          "HALLMARK_ADIPOGENESIS", 
                          "regulation of morphogenesis of an epithelium", 
                          "positive regulation of morphogenesis of an epithelium", 
                          "organic hydroxy compound metabolic process", 
                          "multicellular organismal homeostasis", 
                          "morphogenesis of a branching structure", 
                          "morphogenesis of a branching epithelium", 
                          "embryonic organ development", 
                          "blood vessel morphogenesis", 
                          "blood vessel development", 
                          "extracellular matrix structural constituent", 
                          "external encapsulating structure", 
                          "cornified envelope")

# row_names_df_to_remove %in% combined_results$Description

# Delete uninformative pathways
combined_results_trimmed <- combined_results[!(combined_results$Description %in% row_names_df_to_remove),]

# Calculate the cumulative sums for the number of items in each set
cumulative_sums <- cumsum(table(factor(combined_results_trimmed$Method, levels = c("GSEA GO CC", "GSEA GO MF", "GSEA GO BP", "GSEA"))))

# Create background data with inverted positions
background_data <- data.frame(
  xmin = c(-Inf, -Inf, -Inf, -Inf),  # x-min for all rectangles
  xmax = c(Inf, Inf, Inf, Inf),      # x-max for all rectangles
  ymin = c(0, cumulative_sums[1] + 0.5, cumulative_sums[2] + 0.5, cumulative_sums[3] + 0.5),  # y-min for each rectangle
  ymax = c(cumulative_sums[1] + 0.5, cumulative_sums[2] + 0.5, cumulative_sums[3] + 0.5, length(combined_results_trimmed$Method) + 0.5),  # y-max for each rectangle
  Method = c("GSEA GO CC", "GSEA GO MF", "GSEA GO BP", "GSEA")  # Set names reversed
)

# Define the colors for each method
set_colors <- c("GSEA" = "#FFCCCC", "GSEA GO BP" = "#CCFFCC", "GSEA GO MF" = "#CCCCFF", "GSEA GO CC" = "#CCFCFF")

combined_results_trimmed$Method_num <- match(combined_results_trimmed$Method, c("GSEA GO CC", "GSEA GO MF", "GSEA GO BP", "GSEA"))

# Trim "HALLMARK_" from pathway names
combined_results_trimmed$Description <- str_replace(combined_results_trimmed$Description, "HALLMARK_", "")

# Visualize GSEA
GSEA_sumary <- ggplot(combined_results_trimmed) +
  geom_point(aes(x = NES, y = reorder(Description, Method_num), size = setSize, color = p.adjust)) +  # Ensure 'setSize' and 'p.adjust' are correctly named
  geom_rect(data = background_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Method), 
            inherit.aes = FALSE, alpha = 0.2) +
  scale_fill_manual(values = set_colors) +
  scale_color_gradient(low = "blue", high = "red") +  # Gradient scale for continuous data
  theme_bw() +
  labs(x = "Normalized Enrichment Score (NES)", 
       y = "Pathway", 
       size = "Gene Set Size", 
       color = "Adjusted P-Value")

# Save the plot
ggsave("Graphs_png/Lcp_GSEA_summary.png", plot = GSEA_sumary, height = 10, width = 10)


# Visualization of counts -------------------------------------------------


# Normalize the data
norm_data <- as.data.frame(counts(dds_Lcp_tum, normalized = TRUE))

# Map ENSEMBL ids to gene symbols
norm_data$refs <- mapIds(org.Hs.eg.db, 
                         keys = rownames(norm_data),
                         column = c('SYMBOL'), 
                         keytype = 'ENSEMBL')

# Remove rows with NA values
norm_data <- norm_data[complete.cases(norm_data), ]

# add unique names
# Address duplicate gene names by making them unique
norm_data$refs = make.names(norm_data$refs, unique=TRUE)
rownames(norm_data) <- norm_data$refs

# Transpose the matrix to have genes as columns and samples as rows
norm_data <- norm_data %>% dplyr::select(-refs)
norm_data <- as.data.frame(t(norm_data))

# Extract isotype information
norm_data$isotype <- str_extract(rownames(norm_data), "(IgA|IgG)")

# Filter the dataframe to include only select genes of interest
norm_data <- norm_data[, colnames(norm_data) %in% c("IL5RA", "FCRL4", "RUNX2", "PDCD1", "TNFSF11", "isotype")]

# Reshape data for ggplot boxplot
norm_data_gg <- norm_data %>% 
  pivot_longer(cols = -c("isotype"), names_to = "Var", values_to = "Val")

# Create the boxplot and save as png file
boxplot <- ggplot(norm_data_gg, aes(x = Var, y = Val, fill = isotype)) +
  geom_boxplot() +
  scale_x_discrete(name = "Gene") +
  scale_y_continuous(name = "Normalized Counts") +
  theme_bw()
ggsave(filename = "Graphs_png/Lcp_boxplot_counts.png", plot = boxplot, device = "png", height=5, width=8)

# Extract patient information
norm_data$patient <- str_extract(rownames(norm_data), "(p_\\d\\d|p_\\d)")
norm_data$patient <- str_replace(norm_data$patient, "p_", "p")

# Create dotplots and save them as png files
IL5RA_dotplot <- ggplot(norm_data, aes(x=factor(patient), y=IL5RA, col = isotype)) +
  geom_point(size = 1) +
  scale_x_discrete(name = "Sample") +
  scale_y_continuous(name = "Normalized Counts") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, hjust=1),
        plot.title.position = "plot",        # Move title to middle of the plot
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("IL5RA")
ggsave(filename = "Graphs_png/Lcp_dotplot_counts_IL5RA.png", plot = IL5RA_dotplot, device = "png", height=5, width=5)

FCRL4_dotplot <- ggplot(norm_data, aes(x=factor(patient), y=FCRL4, col = isotype)) +
  geom_point(size = 1) +
  scale_x_discrete(name = "Sample") +
  scale_y_continuous(name = "Normalized Counts") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, hjust=1),
        plot.title.position = "plot",        # Move title to middle of the plot
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("FCRL4") 
ggsave(filename = "Graphs_png/Lcp_dotplot_counts_FCRL4.png", plot = FCRL4_dotplot, device = "png", height=5, width=5)

RUNX2_dotplot <- ggplot(norm_data, aes(x=factor(patient), y=RUNX2, col = isotype)) +
  geom_point(size = 1) +
  scale_x_discrete(name = "Sample") +
  scale_y_continuous(name = "Normalized Counts") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, hjust=1),
        plot.title.position = "plot",        # Move title to middle of the plot
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("RUNX2") 
ggsave(filename = "Graphs_png/Lcp_dotplot_counts_RUNX2.png", plot = RUNX2_dotplot, device = "png", height=5, width=5)

PDCD1_dotplot <- ggplot(norm_data, aes(x=factor(patient), y=PDCD1, col = isotype)) +
  geom_point(size = 1) +
  scale_x_discrete(name = "Sample") +
  scale_y_continuous(name = "Normalized Counts") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, hjust=1),
        plot.title.position = "plot",        # Move title to middle of the plot
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("PDCD1") 
ggsave(filename = "Graphs_png/Lcp_dotplot_counts_PDCD1.png", plot = PDCD1_dotplot, device = "png", height=5, width=5)

RANKL_dotplot <- ggplot(norm_data, aes(x=factor(patient), y=TNFSF11, col = isotype)) +
  geom_point(size = 1) +
  scale_x_discrete(name = "Sample") +
  scale_y_continuous(name = "Normalized Counts") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, hjust=1),
        plot.title.position = "plot",        # Move title to middle of the plot
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("RANKL") 
ggsave(filename = "Graphs_png/Lcp_dotplot_counts_RANKL.png", plot = RANKL_dotplot, device = "png", height=5, width=5)


# # Enrichment of geonomic regions ------------------------------------------
# 
# library(decoupleR)
# library(OmnipathR)
# library(tibble)
# library(dplyr)
# 
# 
# net <- decoupleR::get_dorothea(levels = c('A', 'B', 'C'))
# net
# 
# 
# deg <- res_table_Lcp_tum
# 
# # Access t-values
# # https://support.bioconductor.org/p/120780/
# deg$t_values <- deg$log2FoldChange / deg$lfcSE
# 
# rownames(deg) <- deg$refs
# 
# # Run wmean
# contrast_acts <- run_wmean(mat=deg[, 't_values', drop=FALSE], net=net, .source='source', .target='target',
#                            .mor='mor', times = 100, minsize = 5)
# contrast_acts
# 
# # Filter norm_wmean
# f_contrast_acts <- contrast_acts %>%
#   filter(statistic == 'norm_wmean') %>%
#   mutate(rnk = NA)
# 
# # Filter top TFs in both signs
# n_tfs <- 25
# msk <- f_contrast_acts$score > 0
# f_contrast_acts[msk, 'rnk'] <- rank(-f_contrast_acts[msk, 'score'])
# f_contrast_acts[!msk, 'rnk'] <- rank(-abs(f_contrast_acts[!msk, 'score']))
# tfs <- f_contrast_acts %>%
#   arrange(rnk) %>%
#   head(n_tfs) %>%
#   pull(source)
# f_contrast_acts <- f_contrast_acts %>%
#   filter(source %in% tfs)
# 
# # Plot
# ggplot(f_contrast_acts, aes(x = reorder(source, score), y = score)) + 
#   geom_bar(aes(fill = score), stat = "identity") +
#   scale_fill_gradient2(low = "darkblue", high = "indianred", 
#                        mid = "whitesmoke", midpoint = 0) + 
#   theme_minimal() +
#   theme(axis.title = element_text(face = "bold", size = 12),
#         axis.text.x = 
#           element_text(angle = 45, hjust = 1, size =10, face= "bold"),
#         axis.text.y = element_text(size =10, face= "bold"),
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank()) +
#   xlab("Pathways")
# 
# 
# tf <- 'STAT2'
# 
# df <- net %>%
#   filter(source == tf) %>%
#   arrange(target) %>%
#   mutate(ID = target, color = "3") %>%
#   column_to_rownames('target')
# 
# inter <- sort(intersect(rownames(deg),rownames(df)))
# df <- df[inter, ]
# df[,c('logfc', 't_values', 'p_value')] <- deg[inter, ]
# df <- df %>%
#   mutate(color = if_else(mor > 0 & t_values > 0, '1', color)) %>%
#   mutate(color = if_else(mor > 0 & t_values < 0, '2', color)) %>%
#   mutate(color = if_else(mor < 0 & t_values > 0, '2', color)) %>%
#   mutate(color = if_else(mor < 0 & t_values < 0, '1', color))
# 
# ggplot(df, aes(x = logfc, y = -log10(p_value), color = color, size=abs(mor))) +
#   geom_point() +
#   scale_colour_manual(values = c("red","royalblue3","grey")) +
#   geom_label_repel(aes(label = ID, size=1)) + 
#   theme_minimal() +
#   theme(legend.position = "none") +
#   geom_vline(xintercept = 0, linetype = 'dotted') +
#   geom_hline(yintercept = 0, linetype = 'dotted') +
#   ggtitle(tf)