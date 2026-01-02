library(readr)
library(readxl)
library(scales)
library(cluster)
library(mclust)
library(clusterProfiler)
library(pheatmap)
library(circlize)
library(ggplot2)
library(ggpubr)
library(igraph)
library(here)

# Set Working Directory ----
setwd(here("ASDpipeline"))

sapply(paste0(here(), '/mycode/', list.files(paste0(here(), '/mycode/'))), source)

# Set Figure Save Path
figpath <- 'figs/6_Figures/'
dir.create(figpath) 

#'#######################################################################
# Load Files ----
#'#######################################################################

# Color Scheme
moduleTypecolors = c('Etiological'='#E66100',
                     'Emergent'='#1A85FF',
                     'Other'='#C5B6E0')

# Gene universe 
highvalselgenes <- readRDS('data/highvalselgenes.rds')

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

# logFC 
DE <- read.csv('./data/haney/DEG.csv', header=TRUE)
rownames(DE) <- DE$ensembl_gene_id

# States
states <- readRDS('final/states/maj_states.rds')

#'#######################################################################
# GO Enrichment Analysis ----
#'#######################################################################
# Number of GO terms to include in bar plot
numTerms <- 15
# GO analysis p-value cut off 
pmax <- 0.05

## By Cluster ----
GO_folder <- paste0('GO_byclust/')
dir.create(paste0(figpath, GO_folder))

# Run GO enrichment analysis on each cluster 
GO_results_list <- list()
# Loop through clusters
for (clust in 1:num_clusts){
  cat(clust)
  clust_genes <-  names(clusters[clusters == clust])
  GO_results <- enrichGO(gene = clust_genes, OrgDb = "org.Hs.eg.db",
                         pvalueCutoff = 0.05,
                         keyType = "ENSEMBL", ont = 'BP', pAdjustMethod = 'BH',
                         universe = highvalselgenes)
  if (nrow(GO_results) > 0) {
    GO_results_list[[paste0('C', clust)]] <- GO_results
  }
}
saveRDS(GO_results_list, paste0(figpath, GO_folder, 'GO_ALL_list.rds'))

# Create GO term barplot for each cluster (with same color/p-val scale)
for (clust in 1:num_clusts) {
  cat(clust)
  GO_results <- GO_results_list[[paste0('C', clust)]]
  if (!is.null(GO_results)) {
    if (nrow(GO_results)>0) {
      pdf(paste0(figpath, GO_folder, "GO_ALL0.05_C",clust,".pdf"), width = 7.5, height = 13.5)
      # plot(barplot(GO_results, showCategory = 15))
      p <- barplot(GO_results, showCategory = numTerms)
      p <- p + 
        scale_fill_continuous(limits = c(0, pmax),
                              low='#e06664', high='#327ebb') +
        theme(legend.position = "none")
      plot(p)
      dev.off()
    }
  }
}

# Color bar 
pdf(paste0(figpath, GO_folder, "GO_ALL0.05_colorbar.pdf"), width = 7.5, height = 13.5)
p <- barplot(GO_results_list[[paste0('C', 1)]], showCategory = numTerms)
p <- p + 
  scale_fill_continuous(limits = c(0, pmax),
                        low='#e06664', high='#327ebb')
plot(p)
dev.off()

## Grouped cluster type ----
GO_clusttype_folder <- paste0('GO_grouped/')
dir.create(paste0(figpath, GO_clusttype_folder))

### Grouped Etiological Clusters ----
# Run GO enrichment analysis on grouped etiological genes
GO_etiol_results <- enrichGO(gene = names(clusters[clusters %in% etiol_clusts]), OrgDb = "org.Hs.eg.db",
                             pvalueCutoff = pmax,
                             keyType = "ENSEMBL", ont = 'BP', pAdjustMethod = 'BH',
                             universe = highvalselgenes)
# Create GO term barplot for grouped emergent genes
pdf(paste0(figpath, GO_clusttype_folder, "GO_BP0.05_etiological.pdf"), width = 7.5, height = 13.5)
p <- barplot(GO_etiol_results, showCategory = numTerms)
p <- p + 
  scale_fill_continuous(limits = c(0, pmax),
                        low='#e06664', high='#327ebb') 
plot(p)
dev.off()

### Grouped Emergent Clusters ----
# Run GO enrichment analysis on grouped emergent genes
GO_emerg_results <- enrichGO(gene = names(clusters[clusters %in% emerg_clusts[emerg_clusts <= 19]]), OrgDb = "org.Hs.eg.db",
                             pvalueCutoff = pmax,
                             keyType = "ENSEMBL", ont = 'BP', pAdjustMethod = 'BH',
                             universe = highvalselgenes)
# Create GO term barplot for grouped emergent genes
pdf(paste0(figpath, GO_clusttype_folder, "GO_BP0.05_emergent.pdf"), width = 7.5, height = 13.5)
p <- barplot(GO_emerg_results, showCategory = numTerms)
p <- p + 
  scale_fill_continuous(limits = c(0, pmax),
                        low='#e06664', high='#327ebb') 
plot(p)
dev.off()

# Save grouped GO enrichment analysis results
saveRDS(list('etiological' = GO_etiol_results,
             'emergent' = GO_emerg_results), 
        paste0(figpath, GO_clusttype_folder, 'ASD_GO_BP0.05_clusttype_list.rds'))

# Load grouped GO enrichment analysis results
GO_etiol_results <- readRDS(paste0(figpath, GO_clusttype_folder, 'ASD_GO_BP0.05_clusttype_list.rds'))$etiological
GO_emerg_results <- readRDS(paste0(figpath, GO_clusttype_folder, 'ASD_GO_BP0.05_clusttype_list.rds'))$emergent
# Save original order of GO terms (by adj p-val)
GO_etiol_order <- rownames(data.frame(GO_etiol_results))
GO_emerg_order <- rownames(data.frame(GO_emerg_results))

### Proportion downregulated by GO term ----
# For each grouped GO term, get the proportion of genes (associated with that GO term) that are downregulated (logFC < 0)

# Etiological
etiol_GO15_list <- list()
etiol_counts_df <- data.frame()
etiol_GO15_logFC_df <- data.frame()
for (i in 1:nrow(GO_etiol_results)) {
  GO_row <- data.frame(GO_etiol_results)[i, ]
  GO_id <- GO_row$ID
  GO_descr <- GO_row$Description
  GO_genes <- unlist(strsplit(GO_row$geneID, '/'))
  if (i <= numTerms) {
    for (GO_gene in GO_genes) {
      if (GO_gene %in% names(etiol_GO15_list)) {
        etiol_GO15_list[[GO_gene]] <- paste0(etiol_GO15_list[[GO_gene]], ' // ', GO_descr)
      } else {
        etiol_GO15_list[[GO_gene]] <- GO_descr
      }
    }
  }
  GO_logFCs <- DE$Whole_Cortex_logFC[match(GO_genes, DE$ensembl_gene_id)]
  
  etiol_GO15_logFC_df <- rbind(etiol_GO15_logFC_df, 
                               data.frame('Gene'=GO_genes,
                                          'GO_Num'=rep(i, length(GO_genes)),
                                          'ID'=rep(GO_id, length(GO_genes)),
                                          'Description'=rep(GO_descr, length(GO_genes)),
                                          'logFC'=GO_logFCs))
  
  GO_logFCs <- GO_logFCs[!is.na(GO_logFCs)]
  etiol_counts_df <- rbind(etiol_counts_df, 
                           data.frame('ID'=GO_id,
                                      'Description'=GO_descr,
                                      'NumUpReg'=sum(GO_logFCs > 0),
                                      'NumDownReg'=sum(GO_logFCs < 0),
                                      'PropDownReg'=mean(GO_logFCs < 0)))
  
}
GO_etiol_results <- merge(GO_etiol_results, etiol_counts_df[, c('ID', 'NumUpReg', 'NumDownReg', 'PropDownReg')], by='ID')
# Order in same way as enrichGO
rownames(GO_etiol_results) <- GO_etiol_results$ID
GO_etiol_results <- GO_etiol_results[GO_etiol_order,]
rownames(GO_etiol_results) <- NULL
write.csv(GO_etiol_results, paste0(figpath, GO_clusttype_folder, 'ASD_GO_BP0.05_etiological.csv'))

GO15_etiol_results <- GO_etiol_results[1:numTerms,]

etiol_GO15_genes <- names(etiol_GO15_list)
etiol_GO15_df <- data.frame('Gene'=etiol_GO15_genes,
                            'logFC'=DE$Whole_Cortex_logFC[match(etiol_GO15_genes, DE$ensembl_gene_id)],
                            'GO'=unlist(etiol_GO15_list),
                            'GOcount'=as.numeric(table(unlist(strsplit(data.frame(GO_etiol_results)$geneID[1:numTerms], '/')))[etiol_GO15_genes]))
write.csv(etiol_GO15_df, paste0(figpath, GO_clusttype_folder, 'ASD_logFC_GO15_etiological.csv'))

# Emergent
emerg_GO15_list <- list()
emerg_counts_df <- data.frame()
for (i in 1:nrow(GO_emerg_results)) {
  GO_row <- data.frame(GO_emerg_results)[i, ]
  GO_id <- GO_row$ID
  GO_descr <- GO_row$Description
  GO_genes <- unlist(strsplit(GO_row$geneID, '/'))
  if (i <= numTerms) {
    for (GO_gene in GO_genes) {
      if (GO_gene %in% names(emerg_GO15_list)) {
        emerg_GO15_list[[GO_gene]] <- paste0(emerg_GO15_list[[GO_gene]], ' // ', GO_descr)
      } else {
        emerg_GO15_list[[GO_gene]] <- GO_descr
      }
    }
  }
  GO_logFCs <- DE$Whole_Cortex_logFC[match(GO_genes, DE$ensembl_gene_id)]
  GO_logFCs <- GO_logFCs[!is.na(GO_logFCs)]
  emerg_counts_df <- rbind(emerg_counts_df, 
                           data.frame('ID'=GO_id,
                                      'Description'=GO_descr,
                                      'NumUpReg'=sum(GO_logFCs > 0),
                                      'NumDownReg'=sum(GO_logFCs < 0),
                                      'PropDownReg'=mean(GO_logFCs < 0)))
}
GO_emerg_results <- merge(GO_emerg_results, emerg_counts_df[, c('ID', 'NumUpReg', 'NumDownReg', 'PropDownReg')], by='ID')
# Order in same way as enrichGO
rownames(GO_emerg_results) <- GO_emerg_results$ID
GO_emerg_results <- GO_emerg_results[GO_emerg_order,]
rownames(GO_emerg_results) <- NULL
write.csv(GO_emerg_results, paste0(figpath, GO_clusttype_folder, 'ASD_GO_BP0.05_emergent.csv'))

GO15_emerg_results <- GO_emerg_results[1:numTerms,]

emerg_GO15_genes <- names(emerg_GO15_list)
emerg_GO15_df <- data.frame('Gene'=emerg_GO15_genes,
                            'logFC'=DE$Whole_Cortex_logFC[match(emerg_GO15_genes, DE$ensembl_gene_id)],
                            'GO'=unlist(emerg_GO15_list),
                            'GOcount'=as.numeric(table(unlist(strsplit(data.frame(GO_emerg_results)$geneID[1:numTerms], '/')))[emerg_GO15_genes]))
write.csv(emerg_GO15_df, paste0(figpath, GO_clusttype_folder, 'ASD_logFC_GO15_emergent.csv'))

# Display top 15 GO terms (by adj p-val)
GO_etiol_results <- read.csv(paste0(figpath, GO_clusttype_folder, 'ASD_GO_BP0.05_etiological.csv'), row.names=1)
GO15_etiol_results <- GO_etiol_results[1:numTerms,]

GO_emerg_results <- read.csv(paste0(figpath, GO_clusttype_folder, 'ASD_GO_BP0.05_emergent.csv'), row.names=1)
GO15_emerg_results <- GO_emerg_results[1:numTerms,]

# Shorten description for label
charLim <- 30

etiol_GO_labels <- c()
for (i in 1:numTerms) {
  GOdescr <- as.character(GO15_etiol_results$Description[i])
  if (nchar(GOdescr) > charLim) {
    cat('longer')
    GO_descr_split <- unlist(strsplit(GOdescr, ' '))
    for (j in rev(1:(length(GO_descr_split)-1))) {
      if (nchar(paste0(GO_descr_split[1:j], collapse=' ')) <= charLim) {
        break
      }
    }
    GOlabel <- paste0(paste0(GO_descr_split[1:j], collapse=' '), '\n', paste0(GO_descr_split[j:length(GO_descr_split)], collapse=' '))
  } else {
    GOlabel <- GOdescr
  }
  etiol_GO_labels <- c(etiol_GO_labels, GOlabel)
}
GO15_etiol_results$Label <- etiol_GO_labels

emerg_GO_labels <- c()
for (i in 1:numTerms) {
  GOdescr <- as.character(GO15_emerg_results$Description[i])
  if (nchar(GOdescr) > charLim) {
    cat('longer')
    GO_descr_split <- unlist(strsplit(GOdescr, ' '))
    for (j in rev(1:(length(GO_descr_split)-1))) {
      if (nchar(paste0(GO_descr_split[1:j], collapse=' ')) <= charLim) {
        break
      }
    }
    GOlabel <- paste0(paste0(GO_descr_split[1:j], collapse=' '), '\n', paste0(GO_descr_split[j:length(GO_descr_split)], collapse=' '))
  } else {
    GOlabel <- GOdescr
  }
  emerg_GO_labels <- c(emerg_GO_labels, GOlabel)
}
GO15_emerg_results$Label <- emerg_GO_labels

# Set Description and Label levels so barplot is ordered by propDownReg
GO15_etiol_results$Description <- factor(
  GO15_etiol_results$Description,
  levels = GO15_etiol_results$Description[order(GO15_etiol_results$PropDownReg, decreasing = TRUE)]
)
GO15_etiol_results$Label <- factor(
  GO15_etiol_results$Label,
  levels = GO15_etiol_results$Label[order(GO15_etiol_results$PropDownReg, decreasing = TRUE)]
)

GO15_emerg_results$Description <- factor(
  GO15_emerg_results$Description,
  levels = GO15_emerg_results$Description[order(GO15_emerg_results$PropDownReg, decreasing = TRUE)]
)
GO15_emerg_results$Label <- factor(
  GO15_emerg_results$Label,
  levels = GO15_emerg_results$Label[order(GO15_emerg_results$PropDownReg, decreasing = TRUE)]
)

pdf(paste0(figpath, GO_clusttype_folder, "GO15_BP0.05_etiological_coldownreg_label.pdf"), width = 7.5, height = 13.5)
p <- ggplot(GO15_etiol_results, aes(x=Count, y=Label, fill=PropDownReg)) + 
  geom_col() + 
  scale_fill_continuous(limits = c(0, 1),
                        low='darkorange', high='darkgreen') +
  labs(x = "Count", y = NULL, fill = "Prop Downregulated") +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(size = 12, 
                                    color = "black"))
plot(p)
dev.off() 

pdf(paste0(figpath, GO_clusttype_folder, "GO15_BP0.05_emergent_coldownreg_label.pdf"), width = 7.5, height = 13.5)
p <- ggplot(GO15_emerg_results, aes(x=Count, y=Label, fill=PropDownReg)) + 
  geom_col() + 
  scale_fill_continuous(limits = c(0, 1),
                        low='darkorange', high='darkgreen') +
  labs(x = "Count", y = NULL, fill = "Prop Downregulated")+
  theme_bw() +
  theme(legend.position = "none", 
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(size = 12, 
                                    color = "black"))
plot(p)
dev.off() 

# Try flipping emergent barplot so it's mirrored
pdf(paste0(figpath, GO_clusttype_folder, "GO15_BP0.05_emergent_coldownreg_flipped.pdf"), width = 7.5, height = 13.5)
p <- ggplot(GO15_emerg_results,
            aes(x = -Count, y = Label, fill = PropDownReg)) +
  geom_col() +
  scale_fill_gradient(limits = c(0, 1), low = "darkorange", high = "darkgreen") +
  scale_x_continuous(labels = function(x) abs(x),
                     limits = c(-46, 0),
                     expand = expansion(mult = c(0, 0.02))) +  # show positive labels
  scale_y_discrete(position = "right") +
  labs(x = "Count", y = NULL, fill = "Prop Downregulated") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(size = 12, color = "black"))
print(p)
dev.off() 

# Colorbar
pdf(paste0(figpath, GO_clusttype_folder, "GO15_BP0.05_colorbar_coldownreg.pdf"), width = 7.5, height = 13.5)
p <- ggplot(GO15_emerg_results, aes(x=Count, y=Description, fill=PropDownReg)) + 
  geom_col() + 
  scale_fill_continuous(limits = c(0, 1),
                        low='darkorange', high='darkgreen') +
  labs(x = "Count", y = NULL, fill = "Prop Downregulated")+
  theme_bw() +
  theme(axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(size = 12, 
                                    color = "black"))
plot(p)
dev.off()

#'#######################################################################
# SynGO Enrichment Analyses ----
#'#######################################################################
# Code to generate input for SynGO plotter website to ensure same scale

## By cluster ----

for (clust in c(etiol_clusts, emerg_clusts)) {
  print(clust)
  for (gene in names(clusters[clusters == clust])) {
    cat(gene, ' ')
  }
  cat('\n')
  print('==============================')
}

domain <- 'CC'
# Define your folder path
folder_path <- paste0(figpath, "SynGO_byclust/")  # change this

# Define the prefix you're looking for
prefix <- "SynGO_geneset_analysis__ASD_clust"

# Set cutoff for color/legend
cutoff <- 15.5

for (clust in c(etiol_clusts, emerg_clusts)) {
  # List all files in the folder that start with the prefix
  # Folders have weird date-time related portions in name
  matching_files <- list.files(path = folder_path, pattern = paste0("^", prefix, clust, "__"), full.names = TRUE)
  
  synGO <- read_excel(paste0(matching_files[1], '/syngo_ontologies_with_annotations_matching_user_input.xlsx'))
  cat('synGO', domain, ': clust', clust, '\n')
  synGO_subset <- synGO[!is.na(synGO[,8]) & (synGO[,3] == domain), c(1, 8)]
  for (i in 1:nrow(synGO_subset)) {
    row <- unlist(synGO_subset[i, ])
    val <- -log10(as.numeric(row[2]))
    if (val >= cutoff) {
      val <- cutoff
    }
    cat(row[1], ' ', val)
    cat('\n')
  }
  cat('GO:0098890', ' ', cutoff)
  cat('\n')
  cat('\n')
}

# Define colors and breakpoints
low_color <- "#9EBCDA"
mid_color <- "#756BB1"
high_color <- "#DE2D26"

# Values
min_val <- -3
mid_val <- 0
max_val <- 3

# Create custom palette
my_palette <- colorRampPalette(c(low_color, mid_color, high_color))(100)

# Value to map
val <- 1.2

# Rescale the value to 0–1 relative to min and max
scaled_val <- rescale(val, from = c(min_val, max_val), to = c(0, 1))

# Index into the palette
color_index <- round(scaled_val * (length(my_palette) - 1)) + 1
draw_color <- my_palette[color_index]

print(draw_color)

image(1:100, 1, as.matrix(1:100), col = my_palette, axes = FALSE, xlab = "", ylab = "")

## Grouped cluster type ----

for (gene in names(clusters[clusters %in% etiol_clusts])) {
  cat(gene, ' ')
}
for (gene in names(clusters[clusters %in% emerg_clusts])) {
  cat(gene, ' ')
}

domain <- 'CC'
etiol_synGO <- read_excel(paste0(figpath, 'SynGO_grouped/',
                                 'SynGO_geneset_analysis__ASD_etiol__2025-11-31 12;59/syngo_ontologies_with_annotations_matching_user_input.xlsx'))
etiol_synGO_subset <- etiol_synGO[!is.na(etiol_synGO[,8]) & (etiol_synGO[,3] == domain), c(1, 8)]

emerg_synGO <- read_excel(paste0(figpath, 'SynGO_grouped/',
                                 'SynGO_geneset_analysis__ASD_emerg__2025-11-31 13;00/syngo_ontologies_with_annotations_matching_user_input.xlsx'))
emerg_synGO_subset <- emerg_synGO[!is.na(emerg_synGO[,8]) & (emerg_synGO[,3] == domain), c(1, 8)]

vals <- c()
for (i in 1:nrow(etiol_synGO_subset)) {
  row <- unlist(etiol_synGO_subset[i, ])
  val <- -log10(as.numeric(row[2]))
  vals <- c(vals, val)
  cat(row[1], ' ', val)
  cat('\n')
}
for (i in 1:nrow(emerg_synGO_subset)) {
  row <- unlist(emerg_synGO_subset[i, ])
  val <- -log10(as.numeric(row[2]))
  vals <- c(vals, val)
  cat(row[1], ' ', val)
  cat('\n')
}
max(vals)
min(vals)

cutoff <- 15

for (i in 1:nrow(etiol_synGO_subset)) {
  if (i == 1) {
    cat('GO:0098833', ' ', cutoff)
    cat('\n')
    cat('GO:0098837', ' ', 0.1)
    cat('\n')
  }
  row <- unlist(etiol_synGO_subset[i, ])
  val <- -log10(as.numeric(row[2]))
  if (val >= cutoff) {
    val <- cutoff
  }
  cat(row[1], ' ', val)
  cat('\n')
}
for (i in 1:nrow(emerg_synGO_subset)) {
  if (i == 1) {
    cat('GO:0098833', ' ', cutoff)
    cat('\n')
    cat('GO:0098837', ' ', 0.1)
    cat('\n')
  }
  row <- unlist(emerg_synGO_subset[i, ])
  val <- -log10(as.numeric(row[2]))
  if (val >= cutoff) {
    val <- cutoff
  }
  cat(row[1], ' ', val)
  cat('\n')
}

#'#######################################################################
# logFC Histograms ----
#'#######################################################################
DEdir_figpath <- paste0(figpath, 'logFC/')
dir.create(DEdir_figpath)

## All Genes ----
allgenes <- DE$ensembl_gene_id 
df = data.frame(logFC = DE$Whole_Cortex_logFC[match(allgenes, DE$ensembl_gene_id)])
ggplot(df, aes(x=logFC)) + 
  geom_histogram(alpha = 0.9, position="dodge")+
  xlim(-0.7,0.7)+
  geom_vline(xintercept =0, color = "black", size=0.7) +
  theme_bw()
ggsave(paste0(DEdir_figpath,'/dir_allgenes.pdf'), units = 'in', width=2, height=2)

## 8k Genes ----
df = data.frame(logFC = DE$Whole_Cortex_logFC[match(highvalselgenes, DE$ensembl_gene_id)])
ggplot(df, aes(x=logFC)) + 
  geom_histogram(alpha = 0.9, position="dodge")+
  xlim(-0.7,0.7)+
  geom_vline(xintercept =0, color = "black", size=0.7) +
  theme_bw()
ggsave(paste0(DEdir_figpath,'/dir_8kgenes.pdf'), units = 'in', width=2, height=2)

mean(df$logFC) # -0.03388621
sd(df$logFC) # 0.101457

## Grouped Clusters ----
etiol_genes <- names(clusters[clusters %in% etiol_clusts])
emerg_genes <- names(clusters[clusters %in% emerg_clusts])
group_logFC_df <- data.frame('Gene'=c(etiol_genes, emerg_genes),
                             'logFC'=DE$Whole_Cortex_logFC[match(c(etiol_genes, emerg_genes), DE$ensembl_gene_id)],
                             'clustType'=c(rep('Etiological', length(etiol_genes)),
                                           rep('Emergent', length(emerg_genes))),
                             'Application'=rep('ASD', length(c(etiol_genes, emerg_genes)))
)
write.csv(group_logFC_df, paste0(DEdir_figpath, 'ASD_group_logFC.csv'))
# Create combined plot in ASD pipeline 6_Figures.R file

## By Cluster ----
panel_width <- 1.5
n_panels <- 4

## Etiological Clusters ----
etiol_p=list()
for(clust in sort(etiol_clusts)){
  clustgenes <- names(clusters[clusters == clust])
  
  df = data.frame(logFC = DE$Whole_Cortex_logFC[match(clustgenes, DE$ensembl_gene_id)],
                  state = as.factor(ifelse(unname(states[clustgenes]) == 4, 3, unname(states[clustgenes]))))
  
  etiol_p[[paste0('C', clust)]]=ggplot(df, aes(x=logFC, fill=state)) + 
    geom_histogram(alpha = 0.9, position="dodge")+
    ggtitle(paste0('C', clust))+
    xlim(-0.7,0.7)+
    geom_vline(xintercept =0, color = "black", size=0.7)+
    scale_fill_manual(breaks = c("1","2","3"),
                      labels=c('etiological', 'emergent', 'other'),
                      values = c('#E66100','#1A85FF','#C5B6E0'))+
    theme_bw()
  plot(etiol_p[[paste0('C', clust)]])
  Sys.sleep(1)
}

etiol_plot <- ggarrange(
  plotlist = etiol_p,
  ncol = length(etiol_clusts),
  nrow = 1,
  widths = 6,
  heights = 1,
  common.legend = TRUE,
  legend = 'none'
)
ggsave(paste0(DEdir_figpath,'/dir_etiol_joint.pdf'), 
       etiol_plot,
       width = panel_width * n_panels + 2,
       height = 2, units = "in")

## Emergent Clusters ----
emerg_p=list()
for(clust in sort(emerg_clusts)){
  clustgenes <- names(clusters[clusters == clust])
  
  df = data.frame(logFC = DE$Whole_Cortex_logFC[match(clustgenes, DE$ensembl_gene_id)],
                  state = as.factor(ifelse(unname(states[clustgenes]) == 4, 3, unname(states[clustgenes]))))
  
  emerg_p[[paste0('C', clust)]]=ggplot(df, aes(x=logFC, fill=state)) + 
    geom_histogram(alpha = 0.9, position="dodge")+
    ggtitle(paste0('C', clust))+
    xlim(-0.7,0.7)+
    geom_vline(xintercept =0, color = "black", size=0.7)+
    scale_fill_manual(breaks = c("1","2","3"),
                      labels=c('etiological', 'emergent', 'other'),
                      values = c('#E66100','#1A85FF','#C5B6E0'))+
    theme_bw()
}

emerg_plot <- ggarrange(
  plotlist = emerg_p,
  ncol = length(emerg_clusts),
  nrow = 1,
  widths = 6,
  heights = 1,
  common.legend = TRUE,
  legend = 'none'
)
ggsave(paste0(DEdir_figpath,'/dir_emerg_joint.pdf'), 
       emerg_plot,
       width = panel_width * n_panels + 2,
       height = 2, units = "in")

# Combined plots
c(length(etiol_p), length(emerg_p))
# Same length so no need to fill

comb_p <- c(etiol_p, emerg_p)
comb_plot <- ggarrange(
  plotlist = comb_p,
  ncol = 4,
  nrow = 2,
  widths = 6,
  heights = 1,
  common.legend = TRUE,
  legend = 'none'
)
ggsave(paste0(DEdir_figpath,'/dir_comb_joint.pdf'), 
       comb_plot,
       width = panel_width * 4 + 1,
       height = 4, units = "in")

#'#######################################################################
# Networks ----
#'#######################################################################

# Load network parameter
best.params <- readRDS('final/bestparams.rds')
best.lambda <- best.params$lambda
print(best.lambda)

# Load network
cleanraw <- readRDS(paste0('final/clean_lambda', best.lambda, '.rds'))
graph <- cleanraw$net
set.seed(16)
layout <- layout_with_fr(graph)

graph_figpath <- paste0(figpath, 'graph/')
dir.create(graph_figpath)

## By Cluster ----
# Reload full clusters (including dropped clusters as well)
full_clusters <- readRDS('final/sim_clusters.rds')
table(full_clusters)

# Cluster grouping
full_clust_groups <- readRDS('final/states/clust_groups.rds')
full_etiol_clusts <- full_clust_groups$etiol_clusts
full_emerg_clusts <- full_clust_groups$emerg_clusts
full_other_clusts <- full_clust_groups$other_clusts
full_num_clusts <- max(full_clusters)
full_clust_order <- c(full_etiol_clusts, full_emerg_clusts, full_other_clusts)

# Cluster Color Palette
c('#A50026', '#DD3D2D', '#F67E4B', '#F8945E',
  '#C2E4EF', '#98CAE1', '#82BFD9', '#6EA6CD', '#4A7BB7', '#364B9A',
  '#FDB366',  '#FEDA8B', '#F5E5B6', '#EAECCC')

clust_colorplate <- c(
  # 4 reds/oranges (etiological)
  '#A50026', '#DD3D2D', '#F67E4B', '#F8945E',
  # 6 blues (emergent)                      
  '#364B9A', '#4A7BB7', '#82BFD9', '#C2E4EF', '#6EA6CD', '#98CAE1',
  # 4 yellows (other)
  '#FEDA8B', '#F5E5B6', '#EAECCC', '#FDB366')

clust_colorplate <- c(
  # 4 reds/oranges (etiological)
  '#A50026', '#DD3D2D', '#F67E4B', '#F8945E',
  # 6 blues (emergent) -- 2 dropped
  '#364B9A', '#4A7BB7', '#82BFD9', '#C2E4EF', "#999999", "#808080",
  # 4 yellows (other) -- 1 dropped
  '#FEDA8B', '#F5E5B6', '#EAECCC', "#BFBFBF")

clust_colorref <- rep(NA, full_num_clusts)

clust_colorref[full_etiol_clusts] <- clust_colorplate[1:length(full_etiol_clusts)]
clust_colorref[full_emerg_clusts] <- clust_colorplate[(length(full_etiol_clusts)+1):(length(full_etiol_clusts)+length(full_emerg_clusts))]
clust_colorref[full_other_clusts] <- clust_colorplate[(length(full_etiol_clusts)+length(full_emerg_clusts)+1):full_num_clusts]

# Get color for each gene based on cluster
clust_colors <- clust_colorref[full_clusters]
names(clust_colors) <- names(full_clusters)
# Make sure genes are in same order as graph nodes
clust_colors <- clust_colors[V(graph)$name]

pdf(paste0(graph_figpath, "graph_fr_clust.pdf"), height = 5, width = 5)
par(mar=c(0,0,0,0)+.1)
plot.igraph(graph, vertex.size=V(graph)$size,
            layout=layout,
            vertex.label=NA,
            vertex.color=clust_colors,
            vertex.frame.color=clust_colors,
            edge.width = 0.5,
            edge.color = 'gray')
dev.off()

# Create legend/color bar for clusters
ord <- full_clust_order

# Matrix with numeric values 1..n (needed by pheatmap)
mat <- matrix(seq_along(ord), ncol = 1)
rownames(mat) <- ord

# Colors in the same order, directly from clust_colorref
ordered_clust_colorref <- clust_colorref[ord]

# Colorbar plot
pheatmap(mat,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = ordered_clust_colorref,
         border_color = NA,
         legend = FALSE,
         show_colnames = FALSE,
         filename = paste0(graph_figpath, "graph_clust_colorbar.pdf"),
         width = 0.5, height = 5)

## Colorbar plot with large clusters only ----
ord <- c(etiol_clusts, emerg_clusts, other_clusts)

# Matrix with numeric values 1..n (needed by pheatmap)
mat <- matrix(seq_along(ord), ncol = 1)
rownames(mat) <- ord

# Colors in the same order, directly from clust_colorref
ordered_clust_colorref <- clust_colorref[ord]

# Colorbar plot
pheatmap(mat,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = ordered_clust_colorref,
         border_color = NA,
         legend = FALSE,
         show_colnames = FALSE,
         filename = paste0(graph_figpath, "graph_clust_colorbar_largeclusts.pdf"),
         width = 0.5, height = 5)

## By State ----
# State Color Palette
state_colorref <- c(moduleTypecolors[['Etiological']], moduleTypecolors[['Emergent']], moduleTypecolors[['Other']], moduleTypecolors[['Other']])
# Get color for each gene based on state
state_colors <- state_colorref[states]
names(state_colors) <- names(states)
# Make sure genes are in same order as graph nodes
state_colors <- state_colors[V(graph)$name]

pdf(paste0(graph_figpath, "graph_fr_state.pdf"), height = 5, width = 5)
par(mar=c(0,0,0,0)+.1)
plot.igraph(graph, vertex.size=V(graph)$size,
            layout=layout,
            vertex.label=NA,
            vertex.color=state_colors,
            vertex.frame.color=state_colors,
            edge.width = 0.5,
            edge.color = 'gray')
dev.off()

#'#######################################################################
# Gandal WGCNA Module Comparison ----
#'#######################################################################

clusters_list <- list()
for (i in 1:num_clusts) {
  clusters_list[[i]] <- names(clusters[clusters == i])
}
unlist(lapply(clusters_list, length))

# Load Gandal WGCNA modules
WGCNAmodules <- read_excel('data/SupplementaryTable5.xlsx',sheet = 'Gene_Level', skip=0, col_names=TRUE)
WGCNAmodules$module <- strtoi(sapply(strsplit(sapply(strsplit(WGCNAmodules$WGCNA_module, '_'), '[[', 1), 'M'), '[[', 2))
WGCNAmodules[, c('ensembl_gene_id', 'external_gene_name', 'WGCNA_module', 'module')]
num_modules <- max(WGCNAmodules$module)
WGCNAmodules_list <- list()
for (i in 1:num_modules) {
  WGCNAmodules_list[[i]] <- WGCNAmodules[WGCNAmodules$module == i,]$ensembl_gene_id
}
unlist(lapply(WGCNAmodules_list, length))

# Compare our clusters and the WGCNA modules using the Jaccard Index
WGCNAcomparejacc <- sapply(1:num_modules, function(x){
  sapply(clusters_list,function(y)sum(y%in%WGCNAmodules_list[[x]])/length(unique(c(y, WGCNAmodules_list[[x]]))))
})
colnames(WGCNAcomparejacc) <- paste0('M',1:num_modules)
rownames(WGCNAcomparejacc) <- paste0('C',1:num_clusts)

WGCNAcelltypes <- read_excel('data/SupplementaryTable6.xlsx',sheet = 'GeneModules', skip=1, col_names=TRUE)
WGCNAcelltypes <- as.matrix(WGCNAcelltypes[,2:9])
rownames(WGCNAcelltypes) = paste0('M',1:num_modules)

annotrow <- data.frame('type' = factor(c(rep('Etiological', length(etiol_clusts)),
                                         rep('Emergent', length(emerg_clusts)),
                                         rep('Other', length(other_clusts))), levels=c('Etiological', 'Emergent', 'Other')),
                       'size' = cut(table(clusters)[clust_order],
                                    breaks=c(0,5,10,20,50,100,200,500,1000),
                                    labels=c('<5', '5~10', '10~20', '20~50','50~100','100~200','200~500','>500')))
rownames(annotrow) <- paste0('C', clust_order)
sizecolo <- colorRampPalette(c("white", "darkgreen"))(length(levels(annotrow$size)))
sizecolorrow <- c('<5'=sizecolo[1], '5~10'=sizecolo[2], 
                  '10~20'=sizecolo[3], '20~50'=sizecolo[4],
                  '50~100'=sizecolo[5],'100~200'=sizecolo[6],
                  '200~500'=sizecolo[7],'>500'=sizecolo[8])

annotcol <- data.frame(cellType = as.factor(sapply(rownames(WGCNAcelltypes),function(a)colnames(WGCNAcelltypes)[which.min(WGCNAcelltypes[a,])])))
rownames(annotcol) = colnames(WGCNAcomparejacc)

cellTypecolors <- c('ExNeuron'='yellow','InNeuron'='orange',
                    'Astrocyte'='purple','OPC'='steelblue','Oligo'='pink',
                    'Microglia'='red','Endothelial'='green','Pericyte'='black')
celltypeord <- c(which(annotcol$cellType=='ExNeuron'), which(annotcol$cellType=='InNeuron'),
                 which(annotcol$cellType=='Astrocyte'),which(annotcol$cellType=='OPC'),which(annotcol$cellType=='Oligo'),
                 which(annotcol$cellType=='Microglia'), which(annotcol$cellType=='Endothelial'), which(annotcol$cellType=='Pericyte')
)
for (i in celltypeord) {
  cat(i, ': ')
  minpval <- min(WGCNAcelltypes[i, ])
  cat(colnames(WGCNAcelltypes)[which(WGCNAcelltypes[i, ] == minpval)], '\n')
}
celltypeord <- c(c(11, 12, 14, 17, 19, 31, 3, 5, 26), # ExNeuron
                 which(annotcol$cellType=='InNeuron'),
                 c(2, 29, 32, 22), # Astrocyte
                 which(annotcol$cellType=='OPC'), which(annotcol$cellType=='Oligo'), which(annotcol$cellType=='Microglia'), 
                 c(13, 34, 7, 10, 24), # Endothelial
                 which(annotcol$cellType=='Pericyte')
)
for (i in celltypeord) {
  cat(i, ': ')
  minpval <- min(WGCNAcelltypes[i, ])
  cat(colnames(WGCNAcelltypes)[which(WGCNAcelltypes[i, ] == minpval)], '\n')
}

pheatmap::pheatmap(WGCNAcomparejacc[clust_order, celltypeord],
                   cluster_rows = FALSE, cluster_cols = FALSE,
                   color = colorRampPalette(c("white", "darkblue"))(100),
                   annotation_row = annotrow,
                   annotation_col = annotcol,
                   annotation_colors = list(type=moduleTypecolors,
                                            cellType=cellTypecolors,
                                            size = sizecolorrow
                   ),
                   cellwidth = 30, cellheight = 25, angle_col = "45",
                   width=20,height=15,
                   filename = paste0(figpath,"WGCNA_compareJaccModules.pdf"),
                   fontsize = 20)
