#'##########################################################################
## Second step in finding marker genes and testing marker gene enrichment in clusters
## Kriegstein/Velmeshev SC dataset and Seurat FindMarkers
#'# mayashen@cmu.edu - Sept 2025
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
library(dplyr)
library(reshape2)
library(ggplot2)
library(stringr)
library(pheatmap)
library(here)

# Set Working Directory ----
setwd(here("ADpipeline"))

sapply(paste0(here(), '/mycode/', list.files(paste0(here(), '/mycode/'))), source)

# Set Figure Save Path
figpath <- 'figs/7_MarkerGenes/'
if (!dir.exists(figpath)) {
  dir.create(figpath)
}

dir.create('final/markergenes/')

# Define important variables ----
cellTypes <- c('intNeu', 'excNeu', 'microglia', 'astrocytes', 'opc', 'oligodendrocytes', 'pericytes', 'vascular_cells')
ages <- c("2nd trimester", "3rd trimester", "0-1 years", "1-2 years", "2-4 years", "4-10 years", "10-20 years", "Adult")

ext <- '_pval0.05_logFC1_stable'

# SECTION 1: LOAD FILES, DEFINE VARS ----
#'##########################################################################
#'##########################################################################
#'##                                                                     ###
#'##                 SECTION 1: LOAD FILES, DEFINE VARS                  ###
#'##                                                                     ###
#'##########################################################################
#'##########################################################################
## Load cleaned marker genes ----
# Obtained using Seurat's FindMarkers + post-processing
final_markergenes <- readRDS('final/markergenes/stable_markergenes.rds')

# Gene universe 
univGenes <- readRDS('data/highvalselgenes.rds')

# Clusters
clusters <- readRDS('final/sim_clusters.rds')
table(clusters)

# Cluster grouping
clust_groups <- readRDS('final/states/clust_groups.rds')
etiol_clusts <- clust_groups$etiol_clusts
emerg_clusts <- clust_groups$emerg_clusts
other_clusts <- clust_groups$other_clusts

# Restrict clusters to clusters >= 20
large_clusts <- strtoi(names(table(clusters))[table(clusters) >= 20])
etiol_clusts <- intersect(clust_groups$etiol_clusts, large_clusts)
emerg_clusts <- intersect(clust_groups$emerg_clusts, large_clusts)
other_clusts <- intersect(clust_groups$other_clusts, large_clusts)
clust_order <- intersect(clust_groups$clust_order, large_clusts)

num_clusts <- max(large_clusts)
num_clusts

# SECTION 3: GENERATE CONFUSION MATRICES, ORs, P-VALS ----
#'##########################################################################
#'##########################################################################
#'##                                                                     ###
#'##         SECTION 3: GENERATE CONFUSION MATRICES, ORs, P-VALS         ###
#'##                                                                     ###
#'##########################################################################
#'##########################################################################

# Generate confusion matrices, odds ratios, and corresponding p-values 
# (via Fisher Test) for modules for each cell type

OR_df <- data.frame(OR=numeric(0), pval=numeric(0), module=character(0), celltype=character(0))

for (cellType in cellTypes) {
  print(cellType)
  markerGenes <- intersect(final_markergenes[[cellType]], highvalselgenes)
  
  mx_df = data.frame(Var1=numeric(0), Var2=numeric(0), value=numeric(0), module=character(0))
  fisher_pvals <- c()
  ORs <- c()
  
  for (clust in 1:num_clusts) {
    moduleGenes <- names(which(clusters == clust))
    
    mx <- matrix(0, 2, 2)
    moduleBools <- highvalselgenes %in% moduleGenes
    markerBools <- highvalselgenes %in% markerGenes
    length(moduleBools) == length(markerBools)
    length(markerBools) == length(highvalselgenes)
    
    a <- sum(moduleBools & markerBools)
    b <- sum(moduleBools & !markerBools)
    c <- sum(!moduleBools & markerBools)
    d <- sum(!moduleBools & !markerBools)
    
    dat <- data.frame(
      "markergene_yes" = c(a, c),
      "markergene_no" = c(b, d),
      row.names = c("Module Gene", "Non-Module Gene"),
      stringsAsFactors = FALSE
    )
    colnames(dat) <- c("Marker Gene", "Non-Marker Gene")
    test <- fisher.test(dat, alternative='greater')
    fisher_pvals <- c(fisher_pvals, round(test$p.value, 10))
    ORs <- c(ORs, test$estimate)
    
    if (test$estimate > 50) {
      print(paste0(cellType, ' + ', clust))
      print(dat)
      print(test)
    }
    mx[1, 1] <- a
    mx[1, 2] <- b
    mx[2, 1] <- c
    mx[2, 2] <- d
    
    mx_df_mi <- melt(mx)
    mx_df_mi$module <- rep(paste0('Module ', clust), 4)
    
    mx_df <- rbind(mx_df, mx_df_mi)
  }
  mx_df$module <- factor(mx_df$module, paste0('Module ', 1:num_clusts))
  
  OR_df <- rbind(OR_df, data.frame(OR=ORs,
                                   pval=fisher_pvals,
                                   module=1:num_clusts,
                                   celltype=rep(cellType, num_clusts)))
  
  newLabs <- c()
  for (i in 1:num_clusts) {
    newLabs <- c(newLabs, paste0(levels(mx_df$module)[i], '\n',
                                 'OR: ', round(ORs[i], 5), '\n',
                                 'p-val: ', round(fisher_pvals[i], 5)))
  }
  newLabs <- setNames(newLabs, nm = levels(mx_df$module))
  ggplot(mx_df, aes(x=Var2, y=Var1, fill=value)) +
    geom_tile(color = "white",
              lwd = 1.5,
              linetype = 1) +
    geom_text(aes(label = round(value, 1))) +
    scale_fill_gradient(low = "white", high = "red") +
    coord_fixed() +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank()) +
    facet_wrap(~module, labeller = as_labeller(newLabs), nrow=3) +
    xlab(paste0(str_to_title(cellType), ' Marker Gene: Y or N')) + ylab('Module Gene: Y or N') + labs(title=paste0("Number of Module vs Marker Genes: ", str_to_title(cellType)), fill="Count")
  ggsave(paste0(figpath, cellType, '_markermodule', ext, '.jpg'),
         scale = 1,
         width = 8,
         height = 6,
         units = "in",
         dpi = 800
  )
}
saveRDS(OR_df, paste0('final/markergenes/seuratCleanedMarkersORs', ext, '.rds'))

# SECTION 4: SELECT SIGNIFICANT MODULES ----
#'##########################################################################
#'##########################################################################
#'##                                                                     ###
#'##                SECTION 4: SELECT SIGNIFICANT MODULES                ###
#'##                                                                     ###
#'##########################################################################
#'##########################################################################

OR_df <- readRDS(paste0('final/markergenes/seuratCleanedMarkersORs', ext, '.rds'))

# Drop modules 20-23 from significance consideration (too small)
length(cellTypes) # 8
# 19 clusters/modules
# 8 x 19 = 152 tests
0.05/152

bonfpval <- 0.05
bonfpval_star <- bonfpval/(length(cellTypes) * 19)
OR_df[OR_df$module <= 19 & OR_df$pval <= bonfpval_star,]

# SECTION 5: CREATE HEATMAP ----
#'##########################################################################
#'##########################################################################
#'##                                                                     ###
#'##                      SECTION 5: CREATE HEATMAP                      ###
#'##                                                                     ###
#'##########################################################################
#'##########################################################################

clustNames <- paste0('C', clust_order)
cellTypeNames <- c('InNeuron', 'ExNeuron', 'Microglia', 'Astrocyte', 'OPC', 'Oligo', 'Pericyte', 'Vascular Cell')

OR_mx <- matrix(0, nrow=length(clustNames), ncol=length(cellTypeNames))
rownames(OR_mx) <- clustNames
colnames(OR_mx) <- cellTypeNames

val_method <- 'OR'

for (i in 1:length(clustNames)) {
  for (j in 1:length(cellTypeNames)) {
    if (val_method == 'OR') {
      OR_mx[i, j] <- OR_df[OR_df$module == clust_order[i] & OR_df$celltype == cellTypes[j],]$OR
    } else if (val_method == 'logpval') {
      OR_mx[i, j] <- -log2(OR_df[OR_df$module == clust_order[i] & OR_df$celltype == cellTypes[j],]$pval)
    } else {
      print('invalid value method, val_method must be OR or logpval')
      break
    }
  }
}
# ## UB if there are infinities in matrix ---- 
maxSigOR <- 100
OR_mx[OR_mx > maxSigOR] <- maxSigOR

# SECTION 3: CREATE HEATMAP PLOT ----
#'##########################################################################
#'##########################################################################
#'##                                                                     ###
#'##                    SECTION 3: CREATE HEATMAP PLOT                   ###
#'##                                                                     ###
#'##########################################################################
#'##########################################################################

## Row annotation ----
# Cluster type and size
annotrow <- data.frame('type' = factor(c(rep('Etiological', length(etiol_clusts)),
                                         rep('Emergent', length(emerg_clusts)),
                                         rep('Other', length(other_clusts))), levels=c('Etiological', 'Emergent', 'Other')),
                       'size' = cut(table(clusters)[clust_order],
                                    breaks=c(0,5,10,20,50,100,200,500,1000),
                                    labels=c('<5', '5~10', '10~20', '20~50','50~100','100~200','200~500','>500')))
rownames(annotrow) <- rownames(OR_mx)
sizecolo = colorRampPalette(c("white", "darkgreen"))(length(levels(annotrow$size)))
sizecolorrow = c('<5'=sizecolo[1], '5~10'=sizecolo[2], 
                 '10~20'=sizecolo[3], '20~50'=sizecolo[4],
                 '50~100'=sizecolo[5],'100~200'=sizecolo[6],
                 '200~500'=sizecolo[7],'>500'=sizecolo[8])
## Column annotation ----
# Cell type color, unnecessary because there is only one column per cell type
annotcol <- data.frame('cellType' = cellTypeNames)
rownames(annotcol) <- colnames(OR_mx)

## Annotation colors ----
cellTypecolors = c('InNeuron'='orange', 'ExNeuron'='yellow',
                   'Microglia'='red', 'Astrocyte'='purple',
                   'Oligo'='pink', 'Pericyte'='blue',
                   'Vascular Cell'='green')
moduleTypecolors = c('Etiological'='#E66100',
                     'Emergent'='#1A85FF',
                     'Other'='#C5B6E0')

## Star significant cell type x clusters ----
# Star (*) cells corresponding to marker gene + cell types with significant p-values
test_labels <- matrix('', nrow=length(clustNames), ncol=length(cellTypeNames)) 
for (i in 1:length(clustNames)) {
  for (j in 1:length(cellTypeNames)) {
    if (OR_df[OR_df$module == clust_order[i] & OR_df$celltype == cellTypes[j],]$pval < bonfpval_star) {
      test_labels[i, j] <- '*'
    }
  }
}

# Generate and save heatmap ----
pheatmap::pheatmap(OR_mx,
                   cluster_rows = FALSE, cluster_cols = FALSE,
                   color = colorRampPalette(c("white", "darkblue"))(100),
                   breaks = seq(0, 50, 0.5),
                   annotation_row = annotrow,
                   annotation_colors = list(size=sizecolorrow,
                                            type=moduleTypecolors),
                   display_numbers = test_labels,
                   fontsize_number=25,
                   number_color='red2',
                   cellwidth = 30, cellheight = 25, angle_col = "45",
                   width=12,height=12,
                   filename = paste0(figpath, val_method, ext, '_bonf', bonfpval, '_red2.pdf'),
                   fontsize = 20)
