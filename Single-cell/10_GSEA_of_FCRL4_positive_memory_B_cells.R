library(singleseqgset)
library(heatmap3)

options(ggrepel.max.overlaps = Inf)
set.seed(1)

# set working directory
wd <- "~/B_Memory_master/Master_Diploma_Private/Output/Single-cell/"
setwd(wd)

# Create directory for output if it doesn't exist
if (!dir.exists("Graphs_png")) {
  dir.create("Graphs_png")
}

load("9_Corrected_memory_B_cells.RData")

diff_genes <- FindMarkers(int_2,  ident.1 = 13, min.pct = 0.1, logfc.threshold = -Inf)

any(duplicated(rownames(diff_genes)))

# Let's create a named vector ranked based on the log2 fold change values
lfc_vector <- diff_genes$avg_log2FC
names(lfc_vector) <- rownames(diff_genes)

# We need to sort the log2 fold change values in descending order here
lfc_vector <- sort(lfc_vector, decreasing = TRUE)

hs_hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H")

gsea_results <- GSEA(
  geneList = lfc_vector, # Ordered ranked gene list
  minGSSize = 25, # Minimum gene set size
  maxGSSize = 500, # Maximum gene set set
  pvalueCutoff = 0.05, # p-value cutoff
  eps = 0, # Boundary for calculating the p value
  seed = TRUE, # Set seed to make results reproducible
  pAdjustMethod = "BH", # Benjamini-Hochberg correction
  TERM2GENE = dplyr::select(
    hs_hallmark_sets,
    gs_name,
    gene_symbol
  )
)

head(gsea_results@result)

gsea_result_df <- data.frame(gsea_results@result)

nes_plot_ifn_gamma <- enrichplot::gseaplot(
  gsea_results,
  geneSetID = "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  title = "INTERFERON_GAMMA_RESPONSE",
  color.line = "#0d76ff"
)


nes_plot_ifn_alpha <- enrichplot::gseaplot(
  gsea_results,
  geneSetID = "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  title = "INTERFERON_ALPHA_RESPONSE",
  color.line = "#0d76ff"
)

# Write to .tsv files
write.table(gsea_results@result, 
            file = "Tables/10_GSEA_results_of_13th_cluster.tsv", 
            sep = "\t", 
            row.names = F, 
            quote = FALSE)

# Visualize GSEA output
for (variable in gsea_results@result$ID) {
  tmp <- enrichplot::gseaplot(
    gsea_results,
    geneSetID = variable,
    title = gsub("^HALLMARK_", "", variable),
    color.line = "#0d76ff"
  )
  
  # Modify the png file name with the specific index
  png_file <- paste("Graphs_png/10_", gsub("^HALLMARK_", "", variable), ".png", sep = "")
  png_file_small <- paste("Graphs_png/10_", gsub("^HALLMARK_", "", variable), "_small.png", sep = "")
  
  # Save the plot to the png file using ggsave()
  ggsave(png_file, plot = tmp, height = 14, width = 10)
  ggsave(png_file_small, plot = tmp, height = 6, width = 4)
}

# Trim "HALLMARK_" from pathway names
combined_results_trimmed$Description <- str_replace(combined_results_trimmed$Description, "HALLMARK_", "")

# Visualize GSEA
GSEA_sumary <- ggplot(combined_results_trimmed) +
  geom_point(aes(x = NES, y = Description, size = setSize, color = p.adjust)) +  # Ensure 'setSize' and 'p.adjust' are correctly named
  scale_color_gradient(low = "blue", high = "red") +  # Gradient scale for continuous data
  theme_bw() +
  labs(x = "Normalized Enrichment Score (NES)", 
       y = "Pathway", 
       size = "Gene Set Size", 
       color = "Adjusted P-Value")

# Save the plot
ggsave("Graphs_png/10_FCRL4_plus_GSEA_summary.png", plot = GSEA_sumary, height = 8, width = 10)
ggsave("Graphs_png/10_FCRL4_plus_GSEA_summary_small.png", plot = GSEA_sumary, height = 4, width = 6)
