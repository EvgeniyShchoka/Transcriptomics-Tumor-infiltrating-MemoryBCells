# Comparative transcriptome analysis of B lymphocyte populations infiltrating tumors and expressing IgA and IgG

## Brief description
This repository contains code, tables, and visualizations for a paper soon to be published by Chudakov Lab. The research focuses on investigating the role of tumour-infiltrating memory B cells in Lung Adenocarcinoma (LUAD) and Kidney Renal Clear Cell Carcinoma (KIRC) progression. The findings shed light on the specific subset of FCRL4-expressing memory B cells, indicating exhausted, chronically antigen-stimulated phenotype.

## Repository structure

### Bulk 

Contains scripts for downstream analysis of IgA+ vs IgG+ memory B cells in LUAD and KIRC. The bulk transcriptome libraries obtained in our laboratory will be made available after the release of the associated article. 
The raw data were processed using STAR and featureCounts. GRCh38.p13 genome assembly and GRCh38.109 gene annotation from ensembl.org were used.

[**1_Lcp.R**](/Bulk/R_scripts/1_Lcp.R) 

Deconvolution was performed for identification of contaminated samples. The results were validated by the expression of CD19 and CD20 (MS4A1) marker genes.
<p align="center">
    <img src="https://github.com/EvgeniyShchoka/Transcriptomics-of-IgA-IgG-TIL-B/blob/master/Bulk/Graphs_png/Lcp_heatmap_xCell.png" style="height: 600px;"/>
    <img src="https://github.com/EvgeniyShchoka/Transcriptomics-of-IgA-IgG-TIL-B/blob/master/Bulk/Graphs_png/Lcp_scatterplot_CD19_vs_CD20.png" style="height: 600px;"/>
    
        

### Single-cell

Contains scripts for clusterization and description of tumour-infiltrating memory B cells in LUAD. The single-cell transcriptome data was obtained from Leader et al. [article](https://github.com/effiken/Leader_et_al).

### TCGA

Contains script for Kaplan-Meier curves representing the influence of absolute and normalized value of FCRL4 gene to the progression of LUAD. 

<img src="https://github.com/EvgeniyShchoka/Transcriptomics-of-IgA-IgG-TIL-B/blob/master/TCGA/Graphs_png/surv_plot_normalized_small.png" width=50% height=50%>


