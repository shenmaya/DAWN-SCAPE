# Load Libraries ----
library(cluster)
library(mclust)
library(clusterProfiler)
library(here)

# Set Working Directory ----
setwd(here("ASDpipeline"))

sapply(paste0(here(), '/mycode/', list.files(paste0(here(), '/mycode/'))), source)

# Set Figure Save Path
figpath <- 'figs/4_Clustering/'
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

# Load bulk gene expression data 
control_dat <- readRDS('data/bulk.rds')
control_datreg <- readRDS('data/bulk_reg.rds')

#'#######################################################################
# Cluster Network ----
#'#######################################################################

## Select Clustering Hyperparameters ----
base_seed <- 42
# Set sequence of resol and k to iterate over
resollist <- c(0.0005,0.001,0.005,0.01,0.05,0.1,0.2,0.3)
klist <- seq(4,30,by=1)
reps <- 100

# Create a data frame with all combinations of parameters
param_grid <- expand.grid(resol = resollist, k = klist)
param_grid$id <- as.integer(rownames(param_grid))

folder_cluster <- paste0('pns_clusterrep/')
dir.create(folder_cluster)
bulk <- control_datreg[rownames(cleanraw$ne),]


# Number of cores
num_cores <- detectCores() - 4
num_cores

# Create cluster
cl <- makeCluster(num_cores)

clusterEvalQ(cl, sapply(paste0(here(), '/mycode/', list.files(paste0(here(), '/mycode/'))), source))

# Parallel computation
t0 <- Sys.time()
all_results <- parLapply(
  cl,
  split(param_grid, seq(nrow(param_grid))),
  function(param, cleanraw, bulk, folder, base_seed, reps) {
    run_cluster_net_reps(
      param,
      cleanraw = cleanraw,
      bulk = bulk,
      folder = folder,
      base_seed = base_seed,
      reps = reps,
      savefigs=F
    )
  },
  cleanraw = cleanraw$net,
  bulk = bulk,
  folder = folder_cluster,
  base_seed = base_seed,
  reps = reps
)
stopCluster(cl)
t1 <- Sys.time()
print(t1-t0) 

# Convert the results to a data frame
results_df <- do.call(rbind, all_results)
rownames(results_df) <- 1:nrow(results_df)

dir.create(paste0(folder_cluster, 'metrics/'))
saveRDS(results_df, paste0(folder_cluster, 'metrics/results_', reps, 'reps.rds'))

#'#######################################################################
# Calculate Average Metric ----
#'#######################################################################
# Metric = silhouette score * modularity * stability

## Modularity ----
# Modularity from run_cluster_net_reps

results_df <- readRDS(paste0(folder_cluster, 'metrics/results_', reps, 'reps.rds'))
head(results_df, 12)

## Silhouette Score ----
cor_dat <- cor(t(bulk))
abscor_dat <- abs(cor_dat)

ansraw <- readRDS(paste0('final/', 'pns_lambda', best.lambda,'.rds'))
weightgraph <- ansraw$coef_graphfinal
rownames(weightgraph) = colnames(weightgraph) = ansraw$genesfinal
weightgraph <- weightgraph[rownames(cleanraw$ne), rownames(cleanraw$ne)]

t0 <- Sys.time()
label_silcor_scores <- c()
for (i in 1:nrow(param_grid)) {
  cat(i, '')
  resol <- param_grid[i,]$resol
  k <- param_grid[i,]$k
  
  set.seed(42)
  for (rep in 1:100) {
    labels <- readRDS(paste0('pns_clusterrep/paramsid', i, '_resol', resol, '_k', k,
                             '/labels_resol', resol, '_k', k, '_rep', rep, '.rds'))$labels
    
    silcor_score <- silhouette(labels, 1-abscor_dat)
    label_silcor_scores <- c(label_silcor_scores, mean(silcor_score[, 3]))
  }
}
t1 <- Sys.time()
print(t1-t0)

results_df <- cbind(results_df, data.frame('silcor'=label_silcor_scores))

saveRDS(results_df,  paste0(folder_cluster, 'metrics/run_metrics_df.rds'))

moduas <- c()
silscor <- c()
ks <- c()
resols <- c()
resol_idxs <- c()
for (i in 1:nrow(param_grid)) {
  resol <- param_grid[i,]$resol
  k <- param_grid[i,]$k
  
  param_results <- results_df[results_df$resol == resol & results_df$k == k,]
  if (nrow(param_results) != reps) {
    print(paste0(resol, ' + ', k))
    stop(paste0('not == ', reps))
  }
  moduas <- c(moduas, mean(param_results$modua))
  silscor <- c(silscor, mean(param_results$silcor))
  ks <- c(ks, k)
  resols <- c(resols, resol)
  resol_idxs <- c(resol_idxs, which(resollist == resol))
}

## Jaccard Index (Stability) ---
jaccard_index <- function(cluster1, cluster2) {
  intersect <- sum(cluster1 == cluster2)
  union <- length(cluster1)
  return(intersect / union)
}

tworuns <- combn(reps, 2)
jaccs_list <- list()
jaccs <- c()
t0 <- Sys.time()
for (i in 1:nrow(param_grid)) {
  cat(i, '')
  resol <- param_grid[i,]$resol
  k <- param_grid[i,]$k
  jaccs_i <- c()
  for (j in 1:ncol(tworuns)) {
    run1 <- tworuns[1, j]
    run2 <- tworuns[2, j]
    labels1 <- readRDS(paste0(folder_cluster, 'paramsid', i, '_resol', resol, '_k', k, 
                              '/labels_resol', resol, '_k', k, '_rep', run1, '.rds'))$labels
    labels2 <- readRDS(paste0(folder_cluster, 'paramsid', i, '_resol', resol, '_k', k, 
                              '/labels_resol', resol, '_k', k, '_rep', run2, '.rds'))$labels
    jacc_ij <- jaccard_index(labels1, labels2)
    jaccs_i <- c(jaccs_i, jacc_ij)
  }
  jaccs_list[[paste0('paramsid', i)]] <- jaccs_i
  jaccs <- c(jaccs, mean(jaccs_i))
}
t1 <- Sys.time()
print(t1-t0)

avg_metric_list <- list()
avg_metric_list$moduas <- moduas
avg_metric_list$silscor <- silscor
avg_metric_list$jaccs <- jaccs

saveRDS(avg_metric_list, paste0(folder_cluster, 'metrics/avg_metrics_list.rds'))

avg_metric_list <- readRDS(paste0(folder_cluster, 'metrics/avg_metrics_list.rds'))
moduas <- avg_metric_list$moduas
silscor <- avg_metric_list$silscor
jaccs <- avg_metric_list$jaccs

#'#######################################################################
# Metric Heatmaps ----
#'#######################################################################

modua_mx <- matrix(moduas, nrow=length(klist), ncol=length(resollist), byrow=TRUE)
rownames(modua_mx) <- paste0('k_', klist)
colnames(modua_mx) <- paste0('resol_', resollist)
jpeg(paste0(folder_cluster, 'metrics/metric_modua.jpg'), width = 1500, height = 2000, res = 300, quality = 100)
pheatmap::pheatmap(modua_mx, 
                   cluster_rows=FALSE, cluster_cols=FALSE, 
                   display_numbers = TRUE, number_color = "black")
dev.off()

sil_mx <- matrix(silscor, nrow=length(klist), ncol=length(resollist), byrow=TRUE)
rownames(sil_mx) <- paste0('k_', klist)
colnames(sil_mx) <- paste0('resol_', resollist)
jpeg(paste0(folder_cluster, 'metrics/silh.jpg'), width = 1500, height = 2000, res = 300, quality = 100)
pheatmap::pheatmap(sil_mx, 
                   cluster_rows=FALSE, cluster_cols=FALSE, 
                   display_numbers = TRUE, number_color = "black")
dev.off()

ari_mx <- matrix(aris, nrow=length(klist), ncol=length(resollist), byrow=TRUE)
rownames(ari_mx) <- paste0('k_', klist)
colnames(ari_mx) <- paste0('resol_', resollist)
jpeg(paste0(folder_cluster, 'metrics/ari.jpg'), width = 1500, height = 2000, res = 300, quality = 100)
pheatmap::pheatmap(ari_mx, 
                   cluster_rows=FALSE, cluster_cols=FALSE, 
                   display_numbers = TRUE, number_color = "black")
dev.off()

# What happens if you use the possible ranges instead of the observed ranges, e.g. -1 to +1 instead of 0.29 to 0.73
sils_01_possrange <- (silscor+1)/2
moduas_01_possrange <- (moduas+1)/2
aris_01_possrange <- (aris+1)/2
modua_silh_ari_01_possrange <- sils_01_possrange*moduas_01_possrange*aris_01_possrange
modua_silh_ari_01_possrange_mx <- matrix(modua_silh_ari_01_possrange, nrow=length(klist), ncol=length(resollist), byrow=TRUE)
rownames(modua_silh_ari_01_possrange_mx) <- paste0('k_', klist)
colnames(modua_silh_ari_01_possrange_mx) <- paste0('resol_', resollist)
jpeg(paste0(folder_cluster, 'metrics/metric_ari_modua_silh_01_possrange_possranges.jpg'), width = 1500, height = 2000, res = 300, quality = 100)
pheatmap::pheatmap(modua_silh_ari_01_possrange_mx, 
                   cluster_rows=FALSE, cluster_cols=FALSE, 
                   display_numbers = TRUE, number_color = "black")
dev.off()


sils_01 <- (silscor-min(silscor))/(max(silscor)-min(silscor))
moduas_01 <- (moduas-min(moduas))/(max(moduas)-min(moduas))
aris_01 <- (aris-min(aris))/(max(aris)-min(aris))
modua_silh_ari_01 <- sils_01*moduas_01*aris_01
modua_silh_ari_01_mx <- matrix(modua_silh_ari_01, nrow=length(klist), ncol=length(resollist), byrow=TRUE)
rownames(modua_silh_ari_01_mx) <- paste0('k_', klist)
colnames(modua_silh_ari_01_mx) <- paste0('resol_', resollist)
jpeg(paste0(folder_cluster, 'metrics/metric_ari_modua_silh_01.jpg'), width = 1500, height = 2000, res = 300, quality = 100)
pheatmap::pheatmap(modua_silh_ari_01_mx, 
                   cluster_rows=FALSE, cluster_cols=FALSE, 
                   display_numbers = TRUE, number_color = "black")
dev.off()

overlaid_hists(list('silh'=silscor, 'modu'=moduas, 'ari'=aris), labels=c('Quality (Silhouette Score)', 'Modularity', 'Stability (ARI)'))

overlaid_hists(list('silh'=sils_01, 'modu'=moduas_01, 'ari'=aris_01), labels=c('Quality Norm', 
                                                                               'Modularity Norm',
                                                                               'Stability Norm'))
## Select clustering parameters ----
max_index <- which.max(modua_silh_ari_01_possrange) 
param_grid[which.max(modua_silh_ari_01),] 

best.resol <- param_grid[max_index,]$resol
best.k <- param_grid[max_index,]$k
best.paramid <- param_grid[max_index,]$id

# Save best network/clustering parameters
best.params <- list()
best.params$lambda <- best.lambda
best.params$k <- best.k
best.params$resol <- best.resol
saveRDS(best.params, 'final/bestparams.rds')

best.params <- readRDS('final/bestparams.rds')
best.k <- best.params$k
best.resol <- best.params$resol

#'#######################################################################
# Final Clustering ----
#'#######################################################################

best.params <- readRDS('final/bestparams.rds')
best.lambda <- best.params$lambda
best.k <- best.params$k
best.resol <- best.params$resol
print(c(best.lambda, best.k, best.resol))

# Pre-calculate all gene pairs once
gene_pairs <- t(combn(rownames(cleanraw$ne), 2))  # Transpose for efficient row-wise access

# Initialize similarity matrix
sim_mx <- matrix(0, nrow = nrow(cleanraw$ne), ncol = nrow(cleanraw$ne))
colnames(sim_mx) <- rownames(sim_mx) <- rownames(cleanraw$ne)

t0 <- Sys.time()
# Loop over runs
for (run in 1:reps) {
  cat(run)
  
  labels <- readRDS(paste0(
    'pns_clusterrep/paramsid', best.paramid, 
    '_resol', best.resol, '_k', best.k, 
    '/labels_resol', best.resol, '_k', best.k, '_rep', run, '.rds'
  ))$labels
  
  # Vectorized approach using apply
  same_cluster <- labels[gene_pairs[, 1]] == labels[gene_pairs[, 2]]
  pairs_in_same_cluster <- gene_pairs[same_cluster, , drop = FALSE]
  
  # Update the similarity matrix efficiently
  for (j in 1:nrow(pairs_in_same_cluster)) {
    g1 <- pairs_in_same_cluster[j, 1]
    g2 <- pairs_in_same_cluster[j, 2]
    sim_mx[g1, g2] <- sim_mx[g1, g2] + 1
    sim_mx[g2, g1] <- sim_mx[g2, g1] + 1
  }
}
t1 <- Sys.time()
print(t1 - t0) 

simprop_mx <- sim_mx/reps
saveRDS(simprop_mx, 'final/simprop_mx.rds')

simprop_mx <- readRDS('final/simprop_mx.rds')

sim_hclust <- hclust(as.dist(1-simprop_mx))

# What if some of the cluster sizes are too small? 
# Then let's keep cutting until we get `best.k` good-sized clusters
# We define too small to be 5 or less
k_clust_bool <- F
k <- best.k
while (!k_clust_bool) {
  sim_clusters <- cutree(sim_hclust, k = k)
  if (sum(table(sim_clusters) > 5) >= best.k) {
    k_clust_bool <- T
    break
  }
  k <- k + 1
}
# For ASD we find that all `best.k` (10) clusters have size > 5 

sim_clusters <- cutree(sim_hclust, k = best.k)
table(sim_clusters)
saveRDS(sim_clusters, 'final/sim_clusters_unordered.rds')

sim_clusters <- readRDS('final/sim_clusters_unordered.rds')

# Plot the dendrogram
plot(sim_hclust, main = "hclust from prop sim mx", sub = "", xlab = "", cex = 0.8, labels=FALSE)
# Highlight the k=`best.k` clusters with rectangles
rect.hclust(sim_hclust, k = best.k, border = "red")

# What if we have too big of clusters? 
# For AD, we have clusters of size 300, 289, 216, 202, 200
# From ASD, decide that 250 is the cutoff for too large, aka cut 335 and 271 but not 198

## Split large clusters ----
# Cluster 4 has 335 genes
clust4_genes <- names(which(sim_clusters == 4))
# hclust object for cluster 4  
clust4_sim_hclust <- hclust(as.dist(1-simprop_mx[clust4_genes, clust4_genes]))
# Plot the sub-dendrogram for cluster 4
plot(clust4_sim_hclust, main = "cluster 4 hclust from prop sim mx", sub = "", xlab = "", cex = 0.8, labels=FALSE)
# Visually observing the dendrogram, 3 seems like a good cut
# Highlight the k=3 clusters with rectangles
rect.hclust(clust4_sim_hclust, k = 3, border = "red")
clust4_sim_clusters <- cutree(clust4_sim_hclust, k = 3)
table(clust4_sim_clusters)

# Cluster 1 has 289 genes
clust1_genes <- names(which(sim_clusters == 1))
# hclust object for cluster 1
clust1_sim_hclust <- hclust(as.dist(1-simprop_mx[clust1_genes, clust1_genes]))
# Plot the sub-dendrogram for cluster 1
plot(clust1_sim_hclust, main = "cluster 1 hclust from prop sim mx", sub = "", xlab = "", cex = 0.8, labels=FALSE)
# Visually observing the dendrogram, 3 seems like a good cut
# Highlight the k=3 clusters with rectangles
rect.hclust(clust1_sim_hclust, k = 3, border = "red")
clust1_sim_clusters <- cutree(clust1_sim_hclust, k = 3)
table(clust1_sim_clusters)

# Integrate new (sub-)clusters into main clustering
sim_clusters_split <- c()
for (gene in names(sim_clusters)) {
  if (gene %in% clust4_genes) {
    sim_clusters_split <- c(sim_clusters_split, max(sim_clusters) + clust4_sim_clusters[gene])
  } else if (gene %in% clust1_genes) {
    sim_clusters_split <- c(sim_clusters_split, max(sim_clusters) + length(unique(clust4_sim_clusters)) + clust1_sim_clusters[gene])
  } else {
    sim_clusters_split <- c(sim_clusters_split, sim_clusters[gene])
  }
}
names(sim_clusters_split) <- names(sim_clusters)
table(sim_clusters_split)

## Reorder and remap ----
clust_table <- table(sim_clusters_split)
clust_mapping <- setNames(seq_along(clust_table), names(clust_table)[order(clust_table, decreasing = TRUE)])
renamed_sim_clusters <- clust_mapping[as.character(sim_clusters_split)]
names(renamed_sim_clusters) <- names(sim_clusters_split)

# Verify the new frequency table
table(renamed_sim_clusters)

## Save clusters ----
clusters <- renamed_sim_clusters
table(clusters)
saveRDS(clusters, paste0('final/sim_clusters.rds'))
