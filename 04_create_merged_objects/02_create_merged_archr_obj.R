library(ArchR)
library(parallel)
library(stringr)
library(dplyr)
library(Matrix)
library(readr)
library(BSgenome.Mmusculus.UCSC.mm10)

data_dir <- '/athena/josefowiczlab/scratch/rer4011/projects/MAZ_andrew_data/cellranger_outs/'
archr_dir <- '/athena/josefowiczlab/scratch/rer4011/projects/MAZ_andrew_data/intermediate_results/merged_ArchR'
filtered_data_path <- "/athena/josefowiczlab/scratch/rer4011/projects/MAZ_andrew_data/intermediate_results/per_sample_Archr/"
proj_name <- 'manually_filtered'
sample_folders <- c('MAZ1','MAZ2','MAZ3','MAZ4','MAZ5','MAZ6','MAZ7','MAZ8')

## ################################################################################################
# Merge all manually filtered cells

dfs <- list()
# Loop through each sample to read its respective CSV and process
for (sample in sample_folders) {
  file_path <- file.path(filtered_data_path, paste0(sample,'_ArchR'), "filtered_cells_metadata.csv")
  df <- read.csv(file_path,row.names="X")
  df$orig.ident <- sample 
  df$orig.barcode <- str_split(rownames(df), "#", simplify = TRUE)[, 2]
  dfs[[sample]] <- df
}
# Combine all data frames into one
merged_df <- bind_rows(dfs)
write_csv(merged_df, "/athena/josefowiczlab/scratch/rer4011/projects/tori_atac_data/results/metadata/02_merged_filtered_cells.csv")


# ################################################################################################
# Arrow files and project 

# configure
addArchRThreads(threads = 25) 
addArchRGenome('mm10')
set.seed(42)

# input file paths
inputFiles <- file.path(data_dir, sample_folders, 'outs', 'fragments.tsv.gz')
names(inputFiles) <- sample_folders
barcodes_list <- split(merged_df$orig.barcode, merged_df$orig.ident)

# set output directory
dir.create(archr_dir)
setwd(archr_dir)

#create arrow files
ArrowFiles <- createArrowFiles(
    inputFiles = c(inputFiles),
    sampleNames = names(inputFiles),
    minTSS = 0,
    minFrags = 1,
    maxFrags = 1e+20,
    validBarcodes = barcodes_list,
    addTileMat = TRUE,
    addGeneScoreMat = TRUE,
    excludeChr = c('chrM'),
    verbose = TRUE)

addArchRThreads(threads = 1)

#create archr project
proj <- ArchRProject(
    ArrowFiles = ArrowFiles, 
    outputDirectory = proj_name,
    copyArrows = TRUE
)

# save project
saveArchRProject(ArchRProj = proj, outputDirectory = proj_name, load = FALSE)

# ################################################################################################
# Preprocesing

# IterativeLSI
proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI")

# save project
saveArchRProject(ArchRProj = proj, outputDirectory = proj_name, load = FALSE)

# Peaks
proj <- addClusters(input = proj, reducedDims = "IterativeLSI")
proj <- addGroupCoverages(proj)
proj <- addReproduciblePeakSet(proj)

# Counts
proj <- addPeakMatrix(proj, ceiling=10^9)

# UMAPs
proj <- addUMAP(proj)

# Save
proj <- saveArchRProject(ArchRProj = proj)

# Add doublet score
proj <- addDoubletScores(
  input = proj,
  k = 10,
  knnMethod = "UMAP",
  LSIMethod = 1
)

# Save the project
saveArchRProject(ArchRProj = proj, outputDirectory = proj_name, load = FALSE)

##############################################################################################
#Additional analysis

# Export
dir.create(sprintf("%s/export", proj_name))
write.csv(getReducedDims(proj), sprintf('%s/export/svd.csv', proj_name), quote=FALSE)
write.csv(getCellColData(proj), sprintf('%s/export/cell_metadata.csv', proj_name), quote=FALSE)
write.csv(getEmbedding(proj), sprintf('%s/export/umap.csv', proj_name), quote=FALSE)


# Gene scores
gene.scores <- getMatrixFromProject(proj)
scores <- assays(gene.scores)[['GeneScoreMatrix']]
dir.create(sprintf("%s/export/gene_scores", proj_name))
writeMM(scores, sprintf('%s/export/gene_scores/scores.mtx', proj_name))
write.csv(colnames(scores), sprintf('%s/export/gene_scores/cells.csv', proj_name), quote=FALSE)
write.csv(rowData(gene.scores)$name, sprintf('%s/export/gene_scores/genes.csv', proj_name), quote=FALSE)


# Peak counts
peaks <- getPeakSet(proj)
peak.counts <- getMatrixFromProject(proj, 'PeakMatrix')

# Reorder peaks
# Chromosome order [This mess is necessary since the peaks get sorted by lexicographical order]
chr_order <- sort(seqlevels(peaks))
reordered_features <- list()
for(chr in chr_order)
    reordered_features[[chr]] = peaks[seqnames(peaks) == chr]
reordered_features <- Reduce("c", reordered_features)

# Export counts
dir.create(sprintf("%s/export/peak_counts", proj_name))
counts <- assays(peak.counts)[['PeakMatrix']]
writeMM(counts, sprintf('%s/export/peak_counts/counts.mtx', proj_name))
write.csv(colnames(peak.counts), sprintf('%s/export/peak_counts/cells.csv', proj_name), quote=FALSE)
names(reordered_features) <- sprintf("Peak%d", 1:length(reordered_features))
write.csv(as.data.frame(reordered_features), sprintf('%s/export/peak_counts/peaks.csv', proj_name), quote=FALSE)

# ################################################################################################

# chromVAR scores
proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbp", name = "Motif")
proj <- addBgdPeaks(proj)
proj <- addDeviationsMatrix(
  ArchRProj = proj,
  peakAnnotation = "Motif",
  force = TRUE
)

proj <- saveArchRProject(ArchRProj = proj)

# Motif scores
motifScores <-   getMatrixFromProject(proj, 'MotifMatrix')
dense_matrix <- as.matrix(assays(motifScores)[["z"]])
scores <- as.data.frame(dense_matrix)
write.csv(scores, quote=FALSE,
         sprintf("%s/export/chromvar_motif_scores.csv", proj_name))

# Save the project
saveArchRProject(ArchRProj = proj, outputDirectory = proj_name, load = FALSE)
