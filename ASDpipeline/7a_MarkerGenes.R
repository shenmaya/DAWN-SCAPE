#'##########################################################################
## Data cleaning and first step in finding marker genes
## Kriegstein/Velmeshev SC dataset and Seurat FindMarkers
## Data source: https://www.biorxiv.org/content/10.1101/2022.10.24.513555v1
#'# Followed by 7b_MarkerGenes.R
## mayashen@cmu.edu - Sept 2025
## NOTES: Recommend running on server (not on personal computer) due to size 
## of dataset
#'##########################################################################

# SECTION 0: LOAD LIBRARIES, CHECK DIRECTORY, ETC ----
#'##########################################################################
#'##########################################################################
#'##                                                                     ###
#'##           SECTION 0: LOAD LIBRARIES, CHECK DIRECTORY, ETC           ###
#'##                                                                     ###
#'##########################################################################
#'##########################################################################

## Load libraries ----
library(Matrix)
library(readr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(stringr)
library(here)

# Set Working Directory ----
setwd(here("ASDpipeline"))
getwd()

# SECTION 1: LOAD DATA ----
#'##########################################################################
#'##########################################################################
#'##                                                                     ###
#'##                       SECTION 1: LOAD DATA                          ###
#'##                                                                     ###
#'##########################################################################
#'##########################################################################
# Define important variables ----
given_cellTypes <- c('intNeu', 'excNeu', 'microglia', 'glial_cells', 'pericytes', 'vascular_cells')
ages <- c("2nd trimester", "3rd trimester", "0-1 years", "1-2 years",  "2-4 years", "4-10 years", "10-20 years", "Adult")

t0 <- Sys.time()
# Load full gene expression matrix (all) ----
dataPath <- paste0('data/', 'all', '/')
start_time <- Sys.time()
data <- readMM(paste0(dataPath, 'matrix.mtx'))
end_time <- Sys.time()
print(end_time - start_time)
print(dim(data)) # 20116 709372

scMetaAll <- read.table(paste0(dataPath, 'meta.tsv'),
                        sep="\t", header=1)
nrow(scMetaAll) # 709372

rownames(data) <- read_tsv(paste0(dataPath, 'features.tsv'), col_names=F)$X1 
colnames(data) <- scMetaAll$Cell_ID
print(rownames(data)[1:5]) # NULL
print(colnames(data)[1:5]) # NULL


t1 <- Sys.time()
t1-t0

# Load each cell type's metadata ----
scMetas <- list()
cells <- c()
scClusts <- list()
for (cellType in given_cellTypes) {
  if (cellType != 'astrocytes') {
    dataPath <- paste0('data/', cellType, '/')
    scMeta <- read.table(paste0(dataPath, 'meta.tsv'),
                         sep="\t", header=1)
    print(cellType)
    # print(colnames(scMeta))
    if (cellType %in% c('intNeu', 'vascular_cells')) {
      scMeta$Cell_ID <- scMeta$cellId
    } 
    
    scMetas[[cellType]] <- scMeta
    cells <- c(cells, scMeta$Cell_ID)
    # print(table(scMeta$Seurat_Clusters, useNA='ifany'))
    print(table(scMeta$Cell_Type, useNA='ifany'))
    scClusts[[cellType]] <- strtoi(names(table(scMeta$Seurat_Clusters)))
  }
}
length(cells) # 701256
sum(unlist(lapply(scMetas, nrow))) # 701256
dim(scMetaAll) # 709372
# We seem to be missing some of the cells... investigation in Section 2

# Are there any repeated/overlapping cells across the cell type datasets?
cellTypePairs <- combn(given_cellTypes, 2)
for (i in 1:ncol(cellTypePairs)) {
  cellType1 <- cellTypePairs[,i][1]
  cellType2 <- cellTypePairs[,i][2]
  
  print(paste0(cellType1, ' + ', cellType2, ': ', 
               length(intersect(scMetas[[cellType1]]$Cell_ID, scMetas[[cellType2]]$Cell_ID))))
}
# But at least there seem to be no more repeating cell types?

# Questions:
# 1. Where are the missing cells? 
# Ambiguous or unable to be clustered cells... don't have a good guess for what their cell type is 
# 2. Why are there so many repeated Seurat clusters in the different cell types? Shouldn't each cell type have a unique (or at least close) set of clusters?
# Answer: Seems like it's two sets of clustering, global clustering to get broad cell types and then another clustering within each broad cell type

# Easiest plan: Just drop/ignore the missing cells, might have to do this if I can't figure out what's going on with the Seurat clusters - no way to use that to map cells to a cell type
# SECTION 2: INVESTIGATE EXTRA CELLS ----
#'##########################################################################
#'##########################################################################
#'##                                                                     ###
#'##                 SECTION 2: INVESTIGATE EXTRA CELLS                  ###
#'##                                                                     ###
#'##########################################################################
#'##########################################################################

length(intersect(cells, scMetaAll$Cell_ID)) # 701256 
length(setdiff(scMetaAll$Cell_ID, cells)) # 8116
# There are 8116 cells which are in the all dataset but not anywhere in the CTS datasets...
scMetaAll[scMetaAll$Cell_ID %in% setdiff(scMetaAll$Cell_ID, cells),]
# There doesn't seem to be much rhyme or reason to which cells are missing
# Potentially just unsure what cell type they are but they do have a Seurat cluster... so shouldn't they just be there?
table(scMetaAll[scMetaAll$Cell_ID %in% setdiff(scMetaAll$Cell_ID, cells),]$Lineage)
# AST     ExNeu GLIALPROG       OPC       OUT      VASC 
# 48        95       146        48      5292      2487 
table(scMetaAll[scMetaAll$Cell_ID %in% setdiff(scMetaAll$Cell_ID, cells),]$Seurat_clusters)
# 0    3    4    6   16   18   23   25   26   33   35   38   39 
# 86    5   48   48    1    1    1  146 2487    1 3314 1768  210 
table(scMetaAll[scMetaAll$Cell_ID %in% setdiff(scMetaAll$Cell_ID, cells),]$Dataset)
# Herring     Ramos   Trevino Velmeshev 
#    2156       970        25      4965 

# Investigating extra cells in full dataset
extraCells <- setdiff(scMetaAll$Cell_ID, unique(cells))
extraCellMeans <- colMeans(data[,extraCells]) # none have mean 0
table(scMetaAll[scMetaAll$cell %in% extraCells,]$cluster)
#   0    4   13   20   21   27   28   29   30
# 190   60   14 1369   15 3733 1755 1346  851

# Decide to just continue with cells for which we have cell type labels

# SECTION 4: SELECT CELLS ----
#'##########################################################################
#'##########################################################################
#'##                                                                     ###
#'##                        SECTION 4: SELECT CELLS                      ###
#'##                                                                     ###
#'##########################################################################
#'##########################################################################

# Existing cell types: intNeu, excNeu, microglia, pericytes, vascular cells
# New/sub cell types: astrocytes, oligodendrocytes
cellTypes <- c('intNeu', 'excNeu', 'microglia', 'astrocytes', 'opc', 'oligodendrocytes', 'pericytes', 'vascular_cells')

# Select final set of cells
# Drop glial_progenitors
fullCells <- c()
fullAges <- c()
fullCellTypes <- c()
fullSubCellTypes <- c()
fullPseudoTimes <- c()
for (cellType in cellTypes) {
  print(cellType)
  
  if (cellType %in% given_cellTypes) {
    scMeta <- scMetas[[cellType]]
  } else if (cellType %in% c('astrocytes', 'opc', 'oligodendrocytes')) {
    scMeta <- scMetas[['glial_cells']]
    if (cellType == 'astrocytes') {
      scMeta <- scMeta[scMeta$Cell_Type %in% c('Fibrous_astrocytes', 'Protoplasmic_astrocytes'),]
    } else if (cellType == 'opc') {
      scMeta <- scMeta[scMeta$Cell_Type  == 'OPC',]
    } else if (cellType == 'oligodendrocytes') {
      scMeta <- scMeta[scMeta$Cell_Type  == 'Oligos',]
    }
  } else {
    stop('invalid cellType')
  }
  
  fullCells <- c(fullCells, scMeta$Cell_ID)
  fullAges <- c(fullAges, scMeta$Age_Range)
  fullPseudoTimes <- c(fullPseudoTimes, scMeta$Pseudotime)
  fullCellTypes <- c(fullCellTypes, rep(cellType, nrow(scMeta)))
  
  if (cellType %in% c('intNeu', 'excNeu', 'astrocytes')) {
    fullSubCellTypes <- c(fullSubCellTypes, scMeta$Cell_Type)
  } else {
    fullSubCellTypes <- c(fullSubCellTypes, rep(cellType, nrow(scMeta)))
  }
}

table(fullCellTypes)
table(fullSubCellTypes)

# Create full metadata 
scMetaFinal <- data.frame(cell = fullCells,
                          age = fullAges,
                          pseudotime = fullPseudoTimes,
                          celltype = fullCellTypes,
                          subcelltype = fullSubCellTypes)
rownames(scMetaFinal) <- scMetaFinal$cell

# Create full sc data
scDataFinal <- data[, scMetaFinal$cell]

# Looking at number of cells per group (age x cell type)/group sizes
groupSizeMx <- matrix(0, length(ages), length(cellTypes))
for (c in 1:length(cellTypes)) {
  cellType <- cellTypes[c]
  for (t in 1:length(ages)) {
    age <- ages[t]
    groupSizeMx[length(ages)+1-t, c] <- nrow(scMetaFinal[scMetaFinal$celltype == cellType & scMetaFinal$age == age,])
  }
}
rownames(groupSizeMx) <- rev(ages)
colnames(groupSizeMx) <- cellTypes
groupSizeMx

# Convert and subset genes ----
library(biomaRt)

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

mapping <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype"),
  filters = "external_gene_name",
  values = rownames(scDataFinal),
  mart = mart
)

ADgenes <- readRDS(here('ADpipeline/highvalselgenes.rds'))
ASDgenes <- readRDS(here('ASDpipeline/highvalselgenes.rds'))

length(intersect(mapping$ensembl_gene_id, ADgenes)) # 7504 genes
length(intersect(mapping$ensembl_gene_id, ASDgenes)) # 7501 genes

# Subset down to common genes (alias)
commonGenes <- intersect(mapping$ensembl_gene_id, unique(c(ADgenes, ASDgenes)))
length(commonGenes) # 8022

common_mapping <- mapping[mapping$ensembl_gene_id %in% commonGenes,]
scDataCommon <- scDataFinal[common_mapping$external_gene_name,]

# Change sc expression matrix gene aliases to ensembl IDs
rownames(scDataCommon) <- common_mapping$ensembl_gene_id[match(rownames(scDataCommon), common_mapping$external_gene_name)]

dir.create('data/final/')
saveRDS(mapping, 'data/final/gene_mapping.rds')
saveRDS(scDataCommon, 'data/final/scDataFinal.rds')
saveRDS(scMetaFinal, 'data/final/scMetaFinal.rds')

