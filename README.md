# RNA-seq-Differential-Gene-Expression-Dashboard-
An interactive  dashboard to show the transcriptomic differences between primary and metastatic breast tumours. 
Features 
Data Ingestion: Load raw count matrix and clinical/sample metadata directly from GEO (GSE316391), with support for custom uploads of other count matrices.
Data Validation: Sample-metadata matching checks, missing value detection, library size warnings with helpful error messages 
Pre-processing pipeline: Low-expression gene filtering and DESeq2 median-of-ratios normalization.
DESeq2 Integration: Full differential expression analysis across three contrasts using the gold-standard DESeq2 method.
GO & KEGG Enrichment: Biological Process GO terms and KEGG pathway enrichment on each contrast's DEG list, run immediately after DESeq2 output is generated 
Interactive Visualizations
Volcano plots with gene labeling
Expression heatmaps with hierarchical clustering
MA plots
PCA plots colored by tissue type
Results Export: Download DEG tables as CSV files per contrast 
