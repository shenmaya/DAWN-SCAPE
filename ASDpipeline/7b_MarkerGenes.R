###########################################################################
## Code to find marker genes for Kriegstein/Velmeshev SC dataset using
## Seurat FindMarkers
## mayashen@cmu.edu - Sept 2025
## NOTES: Recommend running on server (not on personal computer) due to size 
## of dataset
###########################################################################

# SECTION 0: LOAD LIBRARIES, CHECK DIRECTORY, ETC ----
#'##########################################################################
#'##########################################################################
#'##                                                                     ###
#'##           SECTION 0: LOAD LIBRARIES, CHECK DIRECTORY, ETC           ###
#'##                                                                     ###
#'##########################################################################
#'##########################################################################

# Load libraries ----
library(SingleCellExperiment)
library(Matrix)
library(plyr)
library(dplyr)
library(Seurat)
library(stringr)
library(here)

# Set Working Directory ----
setwd(here("ASDpipeline"))
getwd()

# Define important variables ----
cellTypes <- c('intNeu', 'excNeu', 'microglia', 'astrocytes', 'opc', 'oligodendrocytes', 'pericytes', 'vascular_cells')
ages <- c("2nd trimester", "3rd trimester", "0-1 years", "1-2 years", "2-4 years", "4-10 years", "10-20 years", "Adult")

# SECTION 1: LOAD DATA ----
#'##########################################################################
#'##########################################################################
#'##                                                                     ###
#'##                         SECTION 1: LOAD DATA                        ###
#'##                                                                     ###
#'##########################################################################
#'##########################################################################

# Load single-cell data ----
scDataFinal <- readRDS('data/final/scDataFinal.rds') # 18135 genes
scMetaFinal <- readRDS('data/final/scMetaFinal.rds')

# SECTION 2: PREPARE VARIABLES ----
#'##########################################################################
#'##########################################################################
#'##                                                                     ###
#'##                     SECTION 2: PREPARE VARIABLES                    ###
#'##                                                                     ###
#'##########################################################################
#'##########################################################################

# Subset cells by age ----
# Want 2 yrs onwards because bulk data is for 2 years onwards
scMetaAge <- scMetaFinal[scMetaFinal$age %in% c('2-4 years', '4-10 years', '10-20 years', 'Adult'),]
scDataAge <- scDataFinal[, rownames(scMetaAge)]
dim(scDataAge)

length(setdiff(colnames(scDataAge), rownames(scMetaAge)))
length(setdiff(rownames(scMetaAge), colnames(scDataAge)))

# Run on ASD and AD genes together as each gene is tested separately
# Less computation time so as to not have to redo the repeated genes


# Create Seurat Object ----
markersSeuratObj<- CreateSeuratObject(counts = scDataAge,
                                      project = "KriegVelmSCE",
                                      meta.data = scMetaAge)

# Set identity class for cell type ----
Idents(object = markersSeuratObj) <- 'celltype'
table(Idents(object = markersSeuratObj))

# Normalize gene expression ----
# markersSeuratObj <- NormalizeData(object = markersSeuratObj, normalization.method = "LogNormalize", scale.factor = 10000)
scDataNorm <- sweep(scDataAge,2,colSums(scDataAge),`/`)
scDataNorm <- log2((scDataNorm*10000) + 1)
scDataNorm <- as(scDataNorm, "CsparseMatrix")

# Add to Seurat object as new assay and set to default assay
markersSeuratObj[["logCPM"]] <- CreateAssayObject(data = scDataNorm)
DefaultAssay(markersSeuratObj) <- "logCPM"
DefaultAssay(markersSeuratObj)

saveRDS(markersSeuratObj, 'data/final/markersSeuratObj.rds')

t0 <- Sys.time()
markersSeuratObj <- readRDS('data/final/markersSeuratObj.rds')
t1 <- Sys.time()
print(t1-t0)

# SECTION 3: FIND MARKER GENES ----
#'##########################################################################
#'##########################################################################
#'##                                                                     ###
#'##                    SECTION 3: FIND MARKER GENES                     ###
#'##                                                                     ###
#'##########################################################################
#'##########################################################################

# Find marker genes using Seurat's FindMarkers
# Each run x cell type gets a new, unique seed 
# Splitting across screens...
# 1:4, seed <- 1
# 5:7, seed <- 33
# 8:10, seed <- 57
seed <- 1
t0 = Sys.time()
for (run in 1:10) {
  print(paste0('run ', run, ', first seed in run: ', seed))
  seuratMarkers = list()
  for (cellType in cellTypes) {
    print(cellType)
    
    geneMarkers <- FindMarkers(markersSeuratObj, ident.1 = cellType, ident.2 = NULL, 
                               only.pos = TRUE, 
                               logfc.threshold=0, 
                               min.diff.pct=0, 
                               max.cells.per.ident=5000,
                               random.seed=seed)
    seuratMarkers[[cellType]] <- geneMarkers
    
    seed <- seed+1
  }
  ext <- paste0('_run', run)
  saveRDS(seuratMarkers, paste0('seuratFullMarkers', ext, '.rds'))
}
t1 = Sys.time()
print(t1 - t0)
print(paste0('next seed: ', seed))

seed <- 1 
t0 = Sys.time()
for (run in 1:10) {
  cat(paste0('run ', run, ', first seed in run: ', seed, '\n'))
  for (cellType in cellTypes) {
    cat(paste0(cellType, '(', seed, ') '))
    seed <- seed+1
  }
  cat('\n =============================================== \n')
}

# SECTION 4: CLEAN SEURAT MARKER GENES ----
#'##########################################################################
#'##########################################################################
#'##                                                                     ###
#'##               SECTION 4: CLEAN SEURAT MARKER GENES                  ###
#'##                                                                     ###
#'##########################################################################
#'##########################################################################
# For each cell type, take potential marker genes with adjusted p-value < 0.05
# For genes that are potential marker genes for multiple cell types, 
# take it to be a potential marker gene for the cell type where it has the largest log2FC
# Then prune potential marker genes by dropping marker genes with log2FC <= 1
# Perform this process 10x with different random seeds (for subsampling), then 
# take marker genes which appear > 5 times in the cell type (majority)

# Back to local 
library(dplyr)
library(here)

# Set Working Directory ----
setwd(here("ASDpipeline"))

cellTypes <- c('intNeu', 'excNeu', 'microglia', 'astrocytes', 'opc', 'oligodendrocytes', 'pericytes', 'vascular_cells')

# Define function to clean markers
getMarkers <- function(markers, pvalcutoff_marker = 0.05, log2FCcutoff_marker=1) {

  ## Select genes with p-value < ?? ----
  markerGenes_df <- data.frame()
  seuratPValMarkers <- list()
  seuratPValMarkers_genes <- list()
  for (cellType in cellTypes) {
    ctMarkers <- markers[[cellType]]
    ctMarkers$ensembl <- rownames(ctMarkers) # FindMarkers
    ctMarkersPVal <- ctMarkers[ctMarkers$p_val_adj < pvalcutoff_marker, ]
    seuratPValMarkers[[cellType]] <- ctMarkersPVal
    seuratPValMarkers_genes[[cellType]] <- rownames(ctMarkersPVal)
    # ctMarkersPVal$ensembl <- rownames(ctMarkersPVal)
    ctMarkersPVal$cellType <- rep(cellType, nrow(ctMarkersPVal))
    markerGenes_df <- rbind(markerGenes_df, ctMarkersPVal)
  }
  
  ## Get unique marker genes by taking max log2FC ----
  markerGenesUnique_df <- markerGenes_df %>%
    # Group the data frame by the unique gene ID
    group_by(ensembl) %>%
    # Within each gene group, select the single row (n = 1) that has the
    # largest (max) value in the avg_log2FC column
    slice_max(order_by = abs(avg_log2FC), n = 1)
  
  ## Drop genes with log2FC <= ?? ----
  markerGenesUniqueFCCutoff_df <- markerGenesUnique_df[abs(markerGenesUnique_df$avg_log2FC) > log2FCcutoff_marker,]
  print(table(markerGenesUniqueFCCutoff_df$cellType))
  
  return(markerGenesUniqueFCCutoff_df)
}

# Collect cleaned potential markers (for each run/seed)
cleaned_markers_list <- list()
for (i in 1:10) {
  print(i)
  markers <- readRDS(paste0('final/markergenes/runs/seuratFullMarkers_run', i, '.rds'))
  markers_cleaned <- getMarkers(markers)
  cleaned_markers_list[[i]] <- markers_cleaned
}

# Collect potential markers within cell types
cellType_markergenes_list <- list()
for (cellType in cellTypes) {
  cat(cellType)
  cellType_markergenes <- c()
  for (i in 1:10) {
    markers_cleaned <- cleaned_markers_list[[i]]
    cellType_markergenes <- c(cellType_markergenes, markers_cleaned[markers_cleaned$cellType == cellType,]$ensembl)
  }
  cellType_markergenes_list[[cellType]] <- cellType_markergenes
  print(table(table(cellType_markergenes)))
}

# Obtain final marker genes by selecting potential markers which appear the majority of the time ( > 5 out of 10 times)
final_markergenes <- list()
for (cellType in cellTypes) {
  cellType_markergenes <- cellType_markergenes_list[[cellType]]
  cellType_final_markergenes <- names(table(cellType_markergenes)[table(cellType_markergenes) > 5])
  cat(paste0(cellType, ': ', length(cellType_final_markergenes)), '\n')
  final_markergenes[[cellType]] <- cellType_final_markergenes
}
saveRDS(final_markergenes, 'final/markergenes/stable_markergenes.rds')
# Save in AD folder as well
saveRDS(final_markergenes, paste0(here('ADpipeline'),
                                  'final/markergenes/stable_markergenes.rds'))
# Note: Need to check respective gene universes (univGenes) within respective applications
