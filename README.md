# Rscripts-ERS
ER stress intersectome in PCOS (analysis code)
This repository contains the R scripts data used to reproduce the main analyses and figures in the manuscript:
The study integrates public GEO datasets (bulk transcriptomics and scRNA-seq) with experimental validation in human granulosa cells to identify and validate hub genes linked to ER stress programs in PCOS.

Contents
R scripts for:
GEO data download and preprocessing (probe annotation, log2 transform, normalization)
Batch correction (ComBat)
Differential expression analysis (limma)
Functional enrichment analysis (GO/KEGG)
Single-gene GSEA
Immune signature scoring (ssGSEA)
Figure generation

Software environment
R version: R software v4.3.3
Recommended OS: Windows/macOS/Linux (tested under standard RStudio workflows)

Key R packages (typical)
The scripts may use (depending on the module):
GEOquery, limma, sva, ggplot2, pheatmap
clusterProfiler, org.Hs.eg.db (or org.Mm.eg.db if needed), enrichplot
GSVA, pROC, rms, rmda (or equivalent DCA package)
Seurat (v5+) and SingleR (for scRNA-seq module)

Public datasets analysed (GEO)
Bulk granulosa cell expression datasets:
GSE34526
GSE10946

Single-cell dataset (human granulosa cells):
GSE240989
All bulk and scRNA-seq datasets are publicly available at NCBI GEO and can be downloaded using GEOquery or via the GEO website.
