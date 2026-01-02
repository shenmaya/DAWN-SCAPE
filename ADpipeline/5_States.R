# Load Libraries ----
library(cluster)
library(mclust)
library(clusterProfiler)
library(dplyr)
library(pheatmap)
library(here)

# Set Working Directory ----
setwd(here("ADpipeline"))

sapply(paste0(here(), '/mycode/', list.files(paste0(here(), '/mycode/'))), source)

# Set Figure Save Path
figpath <- 'figs/5_States/'
if (!dir.exists(figpath)) {
  dir.create(figpath)
}

#'#######################################################################
# Load Files ----
#'#######################################################################

# Load network parameter
best.params <- readRDS('final/bestparams.rds')
best.lambda <- best.params$lambda
print(best.lambda)

# Load network
cleanraw <- readRDS(paste0('final/clean_lambda', best.lambda, '.rds'))
genes <- rownames(cleanraw$ne)
length(genes)
  
# Load final states from runs
final_activegenes <- readRDS(paste0('final/states/final_activegenes.rds'))
final_genes <- readRDS(paste0('final/states/final_genes.rds'))

# Load labels
clusters <- readRDS(paste0('final/sim_clusters.rds'))
table(clusters)

#'#######################################################################
# Get Final Gene States ----
#'#######################################################################

# Setting gene states based on majority vote
maj_states <- c()
for (gene in genes) {
  possible_states <- sapply(final_genes, function(x) x[gene])
  # If tie, defaults to lower state 
  maj_states <- c(maj_states, as.integer(names(which.max(table(possible_states)))))
  
  if (length(unique(possible_states)) > 1) {
    cat(possible_states, ' -> ')
    cat(tail(maj_states, 1))
    cat('\n')
  }
}
names(maj_states) <- genes
table(maj_states)

saveRDS(maj_states, paste0('final/states/maj_states.rds'))

#'#######################################################################
# Create Clust-State Matrix ----
#'#######################################################################

# Value (i, j) is the number of state i genes in cluster j
states <- readRDS(paste0('final/states/maj_states.rds'))

num_active <- c()
num_reactive <- c()
num_other <- c()
num_3s <- c()
num_4s <- c()
num_total <- unname(table(clusters))
for (clust in sort(unique(clusters))){
  cat('Cluster:', clust)
  clust_genes <-  names(clusters[clusters == clust])
  clust_states <- states[clust_genes]
  print(table(clust_states))
  
  num_active <- c(num_active, sum(clust_states==1))
  num_reactive <- c(num_reactive, sum(clust_states==2))
  num_other <- c(num_other, sum(clust_states %in% c(3, 4)))
  num_3s <- c(num_3s, sum(clust_states==3))
  num_4s <- c(num_4s, sum(clust_states==4))
}
clust_state_mx <- matrix(c(num_active, num_reactive, num_other), nrow=3, byrow=TRUE)
rownames(clust_state_mx) <- c('Active', 'Reactive', 'Other')
colnames(clust_state_mx) <- paste0('C', 1:ncol(clust_state_mx))
clust_state_mx
# Save cluster-state matrix 
saveRDS(clust_state_mx, 'final/states/clust_state_mx.rds')

clust_state_mx <- readRDS('final/states/clust_state_mx.rds')

# How to scale cluster-state matrix? Scale by total number of genes within cluster
clust_state_mx_scl_byclust <- sweep(clust_state_mx, 2, colSums(clust_state_mx), '/')
clust_state_mx_scl_byclust

clust_state_mx_scl_bystate <- sweep(clust_state_mx, 1, rowSums(clust_state_mx), '/')
clust_state_mx_scl_bystate

p1 <- pheatmap(t(clust_state_mx_scl_byclust), cluster_rows=F, cluster_cols=F, scale='none', 
               display_numbers=T, color = colorRampPalette(c("white", "red3"))(100), border_color='black', number_color = "black")
p2 <- pheatmap(t(clust_state_mx_scl_bystate), cluster_rows=F, cluster_cols=F, scale='none', 
               display_numbers=T, color = colorRampPalette(c("white", "red3"))(100), border_color='black', number_color = "black")
grid.arrange(p1$gtable, p2$gtable, ncol = 2)
dev.off()

#'#######################################################################
# Classify Clusters ----
#'#######################################################################
etiol_clusts <- c(1, 2, 4, 5, 7, 10, 13, 17)
emerg_clusts <- c(3, 6, 8, 14, 15, 19, 21, 23)
other_clusts <- c(9, 11, 12, 16, 18, 20, 22)
clust_order <- c(other_clusts, emerg_clusts, etiol_clusts)
all(sort(clust_order) == sort(unique(clusters)))
clust_groups <- list()
clust_groups$etiol_clusts <- etiol_clusts
clust_groups$emerg_clusts <- emerg_clusts
clust_groups$other_clusts <- other_clusts
clust_groups$clust_order <- c(etiol_clusts, emerg_clusts, other_clusts)
saveRDS(clust_groups, 'final/states/clust_groups.rds')

clust_groups <- readRDS('final/states/clust_groups.rds')
# Restrict clusters to clusters >= 20
large_clusts <- strtoi(names(table(clusters))[table(clusters) >= 20])
clusters <- clusters[clusters %in% large_clusts]
etiol_clusts <- intersect(clust_groups$etiol_clusts, large_clusts)
emerg_clusts <- intersect(clust_groups$emerg_clusts, large_clusts)
other_clusts <- intersect(clust_groups$other_clusts, large_clusts)
clust_order <- intersect(clust_groups$clust_order, large_clusts)

#'#######################################################################
# Color by Rule ----
#'#######################################################################

rule_dbl_mx <- matrix(0, nrow=nrow(clust_state_mx)*2, ncol=ncol(clust_state_mx))
for (i in 1:23) {
  clust_state_mx_scl_bystate
  clust_state_mx_scl_byclust
  
  if ((clust_state_mx_scl_bystate[1, i] > 0.05) | (clust_state_mx_scl_byclust[1, i] > 0.1)) {
    if (clust_state_mx_scl_bystate[1, i] > 0.05) {
      rule_dbl_mx[1, i] <- 1
    } 
    if (clust_state_mx_scl_byclust[1, i] > 0.1){
      rule_dbl_mx[2, i] <- 1
    }
  } else if ((clust_state_mx_scl_bystate[2, i] > 0.05) | (clust_state_mx_scl_byclust[2, i] > 0.15)) {
    if (clust_state_mx_scl_bystate[2, i] > 0.05) {
      rule_dbl_mx[3, i] <- 1
    } 
    if (clust_state_mx_scl_byclust[2, i] > 0.15) {
      rule_dbl_mx[4, i] <- 1
    }
  } else {
    rule_dbl_mx[5, i] <- 1
    rule_dbl_mx[6, i] <- 1
  }
}
rownames(rule_dbl_mx) <- c('Active - State', 'Active - Clust', 'Reactive - State', 'Reactive - State', 'Other - State', 'Other - Clust')
colnames(rule_dbl_mx) <- paste0('C', 1:ncol(clust_state_mx))

comb_clust_state_mx <- matrix(0, nrow=nrow(clust_state_mx)*2, ncol=ncol(clust_state_mx))
comb_clust_state_mx[1, ] <- clust_state_mx_scl_bystate[1,]
comb_clust_state_mx[3, ] <- clust_state_mx_scl_bystate[2,]
comb_clust_state_mx[5, ] <- clust_state_mx_scl_bystate[3,]
comb_clust_state_mx[2, ] <- clust_state_mx_scl_byclust[1,]
comb_clust_state_mx[4, ] <- clust_state_mx_scl_byclust[2,]
comb_clust_state_mx[6, ] <- clust_state_mx_scl_byclust[3,]

pdf(paste0(figpath, "cluststate_heatmap_dblrule.pdf"), width = 5.625, height = 5.25)
pheatmap(t(rule_dbl_mx[,clust_order]), cluster_rows=F, cluster_cols=F, scale='none', 
         display_numbers=round(t(comb_clust_state_mx[,clust_order]), 2), 
         color = colorRampPalette(c("white", "indianred"))(100), 
         border_color='black', number_color = "black",
         legend=F)
dev.off()

