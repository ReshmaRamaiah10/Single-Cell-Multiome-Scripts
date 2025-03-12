# Multiome (scRNA + scATAC) Analysis Pipeline

This repository contains scripts for analyzing single-cell RNA sequencing (scRNA-seq) and single-cell ATAC sequencing (scATAC-seq) multiome data. The workflow is designed to ensure high-quality cell filtering, annotation, differential analysis, and advanced computational techniques for exploring regulatory mechanisms.

## Table of Contents
1. Installation & Dependencies
2. Data Processing
3. Quality Control & Filtering
4. Cell Annotation
5. Exploratory Data Analysis & Visualization
6. Differential Analysis
7. Advanced Analysis
8. Results Interpretation

### 1. Installation & Dependencies

To run this analysis, install the following software:

CellRanger ARC: Processes raw sequencing data from 10X Genomics into feature-barcode matrices.
Scanpy: Handles scRNA-seq analysis, including normalization, clustering, and differential expression.
ArchR: Processes scATAC-seq data, including peak calling, motif enrichment, and differential accessibility analysis.
Seurat & Signac: Used for multiome integration and additional analysis steps.

### 2. Data Processing

2.1. Running CellRanger ARC
Input raw FASTQ files into CellRanger ARC.
Outputs include a filtered gene expression matrix (for scRNA-seq) and fragment files (for scATAC-seq).
Ensures proper alignment to a reference genome and assigns unique molecular identifiers (UMIs).

### 3. Quality Control & Filtering

3.1. scRNA-seq Quality Control (QC) using Scanpy
Filter out low-quality cells based on:
Number of detected genes per cell.
Total number of UMI counts.
Percentage of mitochondrial gene expression.
Remove doublets using tools like Scrublet or DoubletFinder.
3.2. scATAC-seq Quality Control using ArchR
Assess TSS Enrichment to evaluate chromatin accessibility signal.
Filter low-quality cells based on:
Total number of fragments.
Percentage of reads in peaks.
Doublet identification using computational methods.

### 4. Cell Annotation

4.1. Annotation of scRNA-seq Cells
Cluster cells using graph-based clustering methods (e.g., Leiden or Louvain).
Identify marker genes for each cluster and manually annotate based on known cell types.
Use reference datasets (such as PBMC Atlas) to guide annotations.
4.2. Annotation of scATAC-seq Cells
Compute gene activity scores from chromatin accessibility data.
Transfer labels from scRNA-seq clusters using integration techniques like Harmony or Seurat.
Identify cell-type-specific open chromatin regions using peak calling.

### 5. Exploratory Data Analysis & Visualization

5.1. Visualizing Quality Control Metrics
Generate violin plots and histograms for cell counts, gene counts, and mitochondrial content.
Identify potential batch effects and outliers.
5.2. Clustering and Dimensionality Reduction
Perform PCA, UMAP, and t-SNE to visualize cell heterogeneity.
Assess clustering consistency across different resolutions.
5.3. Cell Type Proportions and Batch Effects
Compare cell-type distribution across different conditions or sample batches.
Apply batch correction techniques if needed (e.g., Harmony, Combat).
5.4. Peak Accessibility Visualization in ATAC Data
Generate genome browser tracks to inspect chromatin accessibility at specific loci.
Plot peak intensities for marker genes associated with different cell types.

### 6. Differential Analysis

6.1. Differential Gene Expression (DGE) Analysis
Compare gene expression between different conditions using Wilcoxon rank-sum test or MAST.
Identify upregulated and downregulated genes across groups.
Perform pseudo-bulk analysis by aggregating cells into sample-level pseudo-replicates.
6.2. Differential Peak Accessibility in ATAC Data
Identify differentially accessible regions (DARs) between conditions.
Perform peak-based differential analysis using statistical tests.
6.3. ChromVAR for Differential Motif Accessibility
Compute TF motif deviations using chromVAR to infer regulatory activity.
Identify transcription factors with altered activity across conditions.

### 7. Advanced Analysis

7.1. Differential Abundance Testing with Milo
Uses MiloR to detect significant changes in the proportion of cell states.
Identifies cell populations expanding or contracting across conditions.
7.2. SPECTRA for Latent Cell State Identification
Uncovers latent regulatory programs using non-negative matrix factorization (NMF).
Identifies transcriptional programs driving cell state changes.
7.3. SCENIC+ for Gene Regulatory Network Inference
Infers gene regulatory networks (GRNs) from multiome data.
Identifies key transcription factors (TFs) and their downstream targets.

### 8. Results Interpretation

8.1. Integrating scRNA-seq and scATAC-seq Insights
Compare gene expression changes with chromatin accessibility patterns.
Identify genes showing coordinated regulation at both the transcript and chromatin level.
8.2. Functional Enrichment Analysis
Perform Gene Ontology (GO) and pathway enrichment analysis for significant genes.
Identify biological processes and pathways associated with differential expression.
8.3. Visualization of Key Findings
Generate summary figures such as heatmaps, volcano plots, and motif enrichment plots.
Create interactive visualizations using tools like cellxgene or UCSC Genome Browser.

### 9. Troubleshooting

Common Issues and Solutions
Issue	Possible Cause	Solution
Low cell yield	Stringent filtering	Adjust QC thresholds carefully
High doublet rate	Poor library quality	Use stricter doublet detection
Batch effects	Technical variability	Apply batch correction techniques
Poor clustering	High noise in data	Optimize dimensionality reduction
