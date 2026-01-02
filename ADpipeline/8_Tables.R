#'##########################################################################
## File containing code to write supplementary tables
#'# mayashen@cmu.edu - Oct 2025
#'# NOTES: 
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

## Load Libraries ----
library(openxlsx)
library(here)

# Set Working Directory ----
setwd(here("ADpipeline"))

sapply(paste0(here(), '/mycode/', list.files(paste0(here(), '/mycode/'))), source)

# Set Figure Save Path
filepath <- 'tables/'
if (!dir.exists(filepath)) {
  dir.create(filepath)
}

filename <- 'AD_supp.xlsx'
# How many tables are there already?
tbl_start <- 10
# How many clusters are there?
clusters <- readRDS('final/sim_clusters.rds')
numAllClusts <- length(table(clusters))
numClusts <- length(table(clusters)[table(clusters) >= 20])
strAllClustsRange <- paste0('(1-', numAllClusts, ')')
strClustsRange <- paste0('(1-', numClusts, ')')

# README Tab ----
# Construct README dataframe as we go
readme_df <- data.frame('TableNumber'=character(),
                        'SheetName'=character(),
                        'Variable'=character(),
                        'Description'=character())

# Table 1: Genes in cleaned network ----
# Load final network
best.params <- readRDS('final/bestparams.rds')
best.lambda <- best.params$lambda
cleanraw <- readRDS(paste0('final/clean_lambda', best.lambda, '.rds'))
genes <- V(cleanraw$net)$name
  
# Load gene information
gen_pvals <- read.csv('data/gen_pvals.csv', row.names=1)
rownames(gen_pvals) <- gen_pvals$ensembl
gene_df <- gen_pvals[genes, c('ensembl', 'alias', 'fixed')]

# Load cluster information 
clusters <- readRDS('final/sim_clusters.rds')
gene_df$cluster <- clusters[genes]

# Load state information
states <- readRDS('final/states/maj_states.rds')
gene_df$state <- states[genes]

colnames(gene_df) <- c('GeneEnsembl', 'GeneSymbol', 'FixedAP', 'Cluster', 'State')

readme_df <- rbind(readme_df, data.frame('TableNumber'=rep(paste0('Table S', tbl_start+1), ncol(gene_df)),
                                         'SheetName'=rep('Genes', ncol(gene_df)),
                                         'Variable'=colnames(gene_df),
                                         'Description'=c("Gene ensembl ID",
                                                         "Gene symbol",
                                                         "Whether the gene is a fixed gene (TRUE=fixed gene, FALSE=not a fixed gene) where the AP p-value is set to 0",
                                                         paste0("Cluster label ", strAllClustsRange),
                                                         "State value")))

# Table 2: Cluster matrix with state type counts ----
# Load cluster matrix with counts of each gene type
clust_mx <- readRDS('final/states/clust_state_mx.rds') 
clust_df <- as.data.frame(t(clust_mx))
clust_df$Cluster <- strtoi(unlist(lapply(strsplit(rownames(clust_df), 'C'), '[[', 2)))
clust_df <- clust_df[, c('Cluster', 'Active', 'Reactive', 'Other')]
colnames(clust_df) <- c('Cluster', 'NumActive', 'NumPReactive', 'NumOther')
clust_df$NumGenes <- rowSums(clust_df[, 2:4])

# Load cluster labels
clust_groups <- readRDS('final/states/clust_groups.rds')
clust_df$ClustLabel <- unlist(lapply(1:numAllClusts, function(x) if (x > numClusts){return(NA)}
                                     else if (x %in% clust_groups$etiol_clusts){return('Etiological')} 
                                     else if (x %in% clust_groups$emerg_clusts){return('Emergent')}
                                     else {'Other'}))

readme_df <- rbind(readme_df, data.frame('TableNumber'=rep(paste0('Table S', tbl_start+2), ncol(clust_df)),
                                         'SheetName'=rep('Clusters', ncol(clust_df)),
                                         'Variable'=colnames(clust_df),
                                         'Description'=c(paste0("Cluster label ", strAllClustsRange),
                                                         "Number of active genes in cluster",
                                                         "Number of p-reactive genes in cluster",
                                                         "Number of other genes in cluster",
                                                         "Number of total genes in cluster",
                                                         "Cluster label (etiological, emergent, or other)")))

# Table 3: Marker genes ----
cellTypes <- c('intNeu', 'excNeu', 'microglia', 'astrocytes', 'opc', 'oligodendrocytes', 'pericytes', 'vascular_cells')
univGenes <- readRDS('data/highvalselgenes.rds')

# Load CTS Marker Genes
markergenes_list <- readRDS(paste0('final/markergenes/stable_markergenes.rds'))
markergenes_df <- data.frame()
for (cellType in cellTypes) {
  ct_markers <- intersect(markergenes_list[[cellType]], univGenes)
  ct_df <- data.frame('CellType'=rep(cellType, length(ct_markers)),
                      'GeneEnsembl'=ct_markers,
                      'GeneAlias'=gen_pvals[ct_markers,]$alias)
  markergenes_df <- rbind(markergenes_df, ct_df)
}

readme_df <- rbind(readme_df, data.frame('TableNumber'=rep(paste0('Table S', tbl_start+3), ncol(markergenes_df)),
                                         'SheetName'=rep('MarkerGenes', ncol(markergenes_df)),
                                         'Variable'=colnames(markergenes_df),
                                         'Description'=c("Cell type for which this is a marker gene for",
                                                         "Gene ensembl ID of marker gene",
                                                         "Gene symbol of marker gene")))

# Table 4: ORs for CTS marker genes and clusters ----
OR_df <- readRDS('final/markergenes/seuratCleanedMarkersORs_pval0.05_logFC1_stable.rds')
colnames(OR_df) <- c('OR', 'PValue', 'Cluster', 'CellType')
OR_df <- OR_df[(OR_df$Cluster <= numClusts),]
OR_df$Significant <- OR_df$PValue < (0.05/nrow(OR_df))
OR_df <- OR_df[, c('OR', 'PValue', 'Significant', 'Cluster', 'CellType')]
readme_df <- rbind(readme_df, data.frame('TableNumber'=rep(paste0('Table S', tbl_start+4), ncol(OR_df)),
                                         'SheetName'=rep('OR', ncol(OR_df)),
                                         'Variable'=colnames(OR_df),
                                         'Description'=c("Odds ratio (OR) for CTS marker gene enrichment in cluster (Fisher's exact test)",
                                                         "P-value for the significance of CTS marker enrichment in cluster (Fisher's exact test)",
                                                         paste0("Cluster label ", strClustsRange),
                                                         "Cell type for which this is a marker gene for",
                                                         "Whether the enrichment for marker genes of the specified cell type is statistically significant (TRUE=significant, FALSE=not significant). Significance is defined by a p-value threshold that has been adjusted for multiple comparisons using the Bonferroni correction (0.05/# tests)")))

# Tables 5+6: GO enrichment by grouped cluster type ----
etiol_df <- read.csv('figs/6_Figures/GO_grouped/AD_GO_BP0.05_etiological.csv', row.names=1)
etiol_df$GeneEnsembl <- etiol_df$geneID
etiol_df$GeneSymbol <- unlist(lapply(etiol_df$geneID, function(x) paste0(gene_df[strsplit(x, '/')[[1]],]$GeneSymbol, collapse='/')))
etiol_df <- etiol_df[, c('ID', 'Description', 'GeneRatio', 'BgRatio', 'RichFactor', 'FoldEnrichment', 'zScore', 'pvalue', 'p.adjust', 'qvalue',
                         'Count', 'NumUpReg', 'NumDownReg', 'PropDownReg',
                         'GeneSymbol', 'GeneEnsembl')]

emerg_df <- read.csv('figs/6_Figures/GO_grouped/AD_GO_BP0.05_emergent.csv', row.names=1)
emerg_df$GeneEnsembl <- emerg_df$geneID
emerg_df$GeneSymbol <- unlist(lapply(emerg_df$geneID, function(x) paste0(gene_df[strsplit(x, '/')[[1]],]$GeneSymbol, collapse='/')))
emerg_df <- emerg_df[, c('ID', 'Description', 'GeneRatio', 'BgRatio', 'RichFactor', 'FoldEnrichment', 'zScore', 'pvalue', 'p.adjust', 'qvalue',
                         'Count', 'NumUpReg', 'NumDownReg', 'PropDownReg',
                         'GeneSymbol', 'GeneEnsembl')]

readme_df <- rbind(readme_df, data.frame('TableNumber'=c(paste0('Table S', tbl_start+7), 
                                                         rep(paste0('Tables S', tbl_start+5, '/S', tbl_start+6, '/S', tbl_start+7), 11), 
                                                         rep(paste0('Tables S', tbl_start+6, '/S', tbl_start+7), 3),
                                                         rep(paste0('Tables S', tbl_start+5, '/S', tbl_start+6, '/S', tbl_start+7), 2)),
                                         'SheetName'=c(rep('GO_Cluster', 1),
                                                       rep('GO_<Etiological/Emergent/Cluster>', 11),
                                                       rep('GO_<Etiological/Emergent>', 3),
                                                       rep('GO_<Etiological/Emergent/Cluster>', 2)),
                                         'Variable'=c('Cluster', colnames(etiol_df)),
                                         'Description'=c(paste0("Cluster label ", strClustsRange),
                                                         "Gene ontology ID",
                                                         "Description of GO term",
                                                         "Number of input genes annotated to GO term / Total number of input genes",
                                                         "Number of background (universe) genes annotated GO term / Total number of background (universe) genes",
                                                         "Number of input genes annotated to GO term / Total number of genes annotated to GO term",
                                                         "GeneRatio/BgRatio",
                                                         "Z-score-like metric representing direction and magnitude of enrichment",
                                                         "Raw p-value from enrichment test",
                                                         "Multiple-testing-corected p-value",
                                                         "Alternative FDR-adjusted p-value",
                                                         "Gene symbols from input list annotated to this GO term, formatted as forward slash-separated list",
                                                         "Gene ensembl IDs from input list annotated to this GO term, formatted as forward slash-separated list",
                                                         "Number of input genes annotated to GO term",
                                                         "Number of up-regulated genes",
                                                         "Number of down-regulated genes",
                                                         "Proportion of down-regulated genes")))

# Table 7: GO enrichment by cluster ----
GOclust_list <- readRDS(paste0('figs/6_Figures/GO_byclust/GO_ALL_list.rds'))

GOclust_df <- data.frame()
for (clusti in 1:numClusts) {
  clust <- paste0('C', clusti)
  if (clust %in% names(GOclust_list)) {
    GOclusti_df <- data.frame(GOclust_list[[clust]])
    GOclusti_df$Cluster <- rep(clusti, nrow(GOclusti_df))
    GOclusti_df$GeneEnsembl <- GOclusti_df$geneID
    GOclusti_df$GeneSymbol <- unlist(lapply(GOclusti_df$geneID, function(x) paste0(gene_df[strsplit(x, '/')[[1]],]$GeneSymbol, collapse='/')))
    GOclusti_df <- GOclusti_df[, c('Cluster', 'ID', 'Description', 'GeneRatio', 'BgRatio', 'RichFactor', 'FoldEnrichment', 'zScore', 'pvalue', 'p.adjust', 'qvalue',
                                   'Count',
                                   'GeneSymbol', 'GeneEnsembl')]
    GOclust_df <- rbind(GOclust_df, GOclusti_df)
  }
}

# Create xlsx object
wb <- openxlsx::createWorkbook()
# Add README sheet
openxlsx::addWorksheet(wb, 'README')
openxlsx::writeData(
  wb, 
  sheet = 'README', 
  x = readme_df, 
  colNames = TRUE, 
  rowNames = FALSE # Matches your original intent
)
# Add Genes sheet
openxlsx::addWorksheet(wb, 'Genes')
openxlsx::writeData(
  wb, 
  sheet = 'Genes', 
  x = gene_df, 
  colNames = TRUE, 
  rowNames = FALSE # Matches your original intent
)
# Add Clusters sheet
openxlsx::addWorksheet(wb, 'Clusters')
openxlsx::writeData(
  wb, 
  sheet = 'Clusters', 
  x = clust_df, 
  colNames = TRUE, 
  rowNames = FALSE # Matches your original intent
)
# Add MarkerGenes sheet
openxlsx::addWorksheet(wb, 'MarkerGenes')
openxlsx::writeData(
  wb, 
  sheet = 'MarkerGenes', 
  x = markergenes_df, 
  colNames = TRUE, 
  rowNames = FALSE # Matches your original intent
)
# Add OR sheet
openxlsx::addWorksheet(wb, 'OR')
openxlsx::writeData(
  wb, 
  sheet = 'OR', 
  x = OR_df, 
  colNames = TRUE, 
  rowNames = FALSE # Matches your original intent
)
# Add GO_Etiological sheet
openxlsx::addWorksheet(wb, 'GO_Etiological')
openxlsx::writeData(
  wb, 
  sheet = 'GO_Etiological', 
  x = etiol_df, 
  colNames = TRUE, 
  rowNames = FALSE # Matches your original intent
)
# Add GO_Emergent sheet
openxlsx::addWorksheet(wb, 'GO_Emergent')
openxlsx::writeData(
  wb, 
  sheet = 'GO_Emergent', 
  x = emerg_df, 
  colNames = TRUE, 
  rowNames = FALSE # Matches your original intent
)
# Add GO_Cluster sheet
openxlsx::addWorksheet(wb, 'GO_Cluster')
openxlsx::writeData(
  wb, 
  sheet = 'GO_Cluster', 
  x = GOclust_df, 
  colNames = TRUE, 
  rowNames = FALSE # Matches your original intent
)
openxlsx::saveWorkbook(wb, paste0(filepath, filename), overwrite = TRUE)

