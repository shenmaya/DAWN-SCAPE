# Load Libraries ----
library(doBy)
library(ggplot2)
library(here)

# Set Working Directory ----
setwd(here("ASDpipeline"))

sapply(paste0(here(), '/mycode/', list.files(paste0(here(), '/mycode/'))), source)

# Set Figure Save Path
figpath <- 'figs/2_ConstructNetwork/'
if (!dir.exists(figpath)) {
  dir.create(figpath)
}

metrics <- c('meanDegree', 'medianDegree', 'R2', 
             'propDeg1', 'meanShortestPath', 'normMeanShortestPath', 
             'meanCVMSE', 'medCVMSE', 'propStblEdges')
metric_titles <- c('Average Degree', 'Median Degree', 'R2', 
                   'Proportion Degree 1', 'Average Shortest Path', 'Normalized Average Shortest Path',
                   'Average Cross-Validated MSE', 'Median Cross-Validated MSE', 'Proportion Stable Edges')

#'#######################################################################
# Load Data ----
#'#######################################################################

# Load z-scores 
gen_df <- read.csv('data/gen_pvals.csv', row.names=1)
de_df <- read.csv('data/de_pvals.csv', row.names=1)
selgenes <- gen_df$ensembl

# Load bulk data
control_dat <- readRDS('data/bulk.rds')
control_datreg <- readRDS('data/bulk_reg.rds')

#'#######################################################################
# Subset Genes ----
#'#######################################################################

# Select genes which are over 50% nonzero (over cells)
highnonzeroselgenes = selgenes[which(rowMeans(control_dat[selgenes,]>0)>0.5)]
str(highnonzeroselgenes) #13522
# From previous set, select 8000 genes with largest mean values
highvalselgenes = highnonzeroselgenes[which.maxn(rowMeans(control_dat[highnonzeroselgenes,]),8000)]
str(highvalselgenes) # 8000
saveRDS(highvalselgenes, 'data/highvalselgenes.rds')

# From previous set, select genes with either TADA p-value or DE p-value below 10th percentile (respectively) or both 
highsubselgenes <- highvalselgenes[which(gen_df$pvals[match(highvalselgenes, gen_df$ensembl)] < quantile(gen_df$pvals[match(highvalselgenes, gen_df$ensembl)],0.1) | 
                                           de_df$pvals[match(highvalselgenes, de_df$ensembl)] < quantile(de_df$pvals[match(highvalselgenes, de_df$ensembl)],0.1))]
str(highsubselgenes) # 1532

#'#######################################################################
# Estimate Network ----
#'#######################################################################

X <- t(control_datreg[highvalselgenes, ])
pvals <- gen_df[match(highvalselgenes, gen_df$ensembl),]$pval
names(pvals) <- highvalselgenes

# For running on server
dir.create('2_ConstructNetwork/')
save(X, pvals, highsubselgenes, file = "2_ConstructNetwork/server_constructnetwork.RData") 

# Server code
codepath <- '/raid6/home/myshen/gene-community-server/mycode/'
sapply(paste0(codepath, list.files(codepath)), source)
server_path <- '2_ConstructNetwork/'
load(file = paste0(server_path, "server_constructnetwork.RData"))

savepath_rd1 <- paste0(server_path, 'stability_rd1/')
dir.create(savepath_rd1)
lambdas_rd1 <-  seq(0.1, 0.6, 0.05)
t0 <- Sys.time()
lambda_stability_graphs_rd1 <- construct_graphs_stability(X, pvals, pvalthresh = NULL, corthresh = 0.5, pvalscreen = highsubselgenes, 
                                                          lambdas = lambdas_rd1, numsmps=100, propsubsmp=0.8,
                                                          degcut=1, maxnumcomp=NULL, mincompsize=10, 
                                                          num_cores=50, savefile=savepath_rd1, standardize=T, verbose=T)
t1 <- Sys.time()
print(t1-t0) # Time difference of 1.590162 hours

savepath_rd2 <- paste0(server_path, 'stability_rd2/')
dir.create(savepath_rd2)
lambdas_rd2 <-  seq(0.15, 0.3, 0.01)
t0 <- Sys.time()
lambda_stability_graphs_rd2 <- construct_graphs_stability(X, pvals, pvalthresh = NULL, corthresh = 0.5, pvalscreen = highsubselgenes, 
                                                          lambdas = lambdas_rd2, numsmps=100, propsubsmp=0.8,
                                                          degcut=1, maxnumcomp=NULL, mincompsize=10, 
                                                          num_cores=50, savefile=savepath_rd2, standardize=T, verbose=T)
t1 <- Sys.time()
print(t1-t0) # Time difference of 1.778603 hours

## Round 1 ----
folder_rd1 <- 'pns_metrics_rd1/'
dir.create(folder_rd1)
vizfolder_rd1 <- paste0(folder_rd1, 'viz/')
dir.create(vizfolder_rd1)
stblfolder_rd1 <- paste0(folder_rd1, 'stability/')
dir.create(stblfolder_rd1)

lambdas_rd1 <-  seq(0.1, 0.6, 0.05)

### Easier Metrics ----
metrics_rd1 <- construct_graphs_metrics(X, pvals, pvalthresh=NULL, corthresh=0.5, 
                                        pvalscreen=highsubselgenes, lambdas=lambdas_rd1, 
                                        degcut=1, maxnumcomp=NULL, mincompsize=10, degdivide=5,
                                        standardize=T, savefile=folder_rd1, saveviz=vizfolder_rd1, verbose=T)
saveRDS(metrics_rd1, paste0(folder_rd1, 'metrics_rd1.rds'))
metrics_rd1 <- readRDS(paste0(folder_rd1, 'metrics_rd1.rds'))

metrics_rd1$pns_normMeanShortestPath <- metrics_rd1[[paste0('pns_', 'meanShortestPath')]]/metrics_rd1[[paste0('pns_', 'numNodes')]]
metrics_rd1$clean_normMeanShortestPath <- metrics_rd1[[paste0('clean_', 'meanShortestPath')]]/metrics_rd1[[paste0('clean_', 'numNodes')]]

### CV Mean Squared Error ----
# Resulting matrix: (# core genes) x (# lambdas)
# PNS: For lambda, for each core gene, run linear regression on other core genes + nbr genes --> error
# MSE_CV: For each lambda, for each core gene, PNS with K-fold CV --> K-Fold CV error for each gene, average over K folds 
# --> Return matrix of mean K-fold CVs, one for each core gene and one for each lambda
t0 <- Sys.time()
cvmse_mx_rd1 <- get_MSEs_CV(X, pvals, pvalthresh=NULL, corthresh=0.5, pvalscreen=highsubselgenes, lambdas=lambdas_rd1, 
                           standardize=T, verbose=T, K=5)
t1 <- Sys.time()
print(t1-t0)
saveRDS(cvmse_mx_rd1, paste0(folder_rd1, 'cvmse_mx_rd1.rds'))

# Save mean and median CV MSE as both pns and clean for easier plotting
mean_cvmse_rd1 <- colMeans(cvmse_mx_rd1)
metrics_rd1$pns_meanCVMSE <- mean_cvmse_rd1
metrics_rd1$clean_meanCVMSE <- mean_cvmse_rd1

med_cvmse_rd1 <- apply(cvmse_mx_rd1, 2, median)
metrics_rd1$pns_medCVMSE <- med_cvmse_rd1
metrics_rd1$clean_medCVMSE <- med_cvmse_rd1

### Stability ----
# Can run on server
numsmps <- 100
prop_thresh <- 0.6
prop_base_thresh <- 0.1
pns_num_stable_edges <- c()
pns_num_conn_nodes <- c()
pns_num_overall_edges <- c()
clean_num_stable_edges <- c()
clean_num_conn_nodes <- c()
clean_num_overall_edges <- c()
for (i in 1:length(lambdas_rd1)) {
  lambda <- lambdas_rd1[i]
  cat(lambda, '')
  pns_agg_mx <- readRDS(paste0(stblfolder_rd1, 'pns_lambda', lambda, '_numsmps', numsmps, '_propsubsmp0.8.rds'))
  pns_agg_prop_mx <- pns_agg_mx/numsmps
  pns_num_stable_edges <- c(pns_num_stable_edges, sum(pns_agg_prop_mx[upper.tri(pns_agg_prop_mx)] > prop_thresh))
  pns_num_overall_edges <- c(pns_num_overall_edges, sum(pns_agg_prop_mx[upper.tri(pns_agg_prop_mx)] > prop_base_thresh))
  pns_num_conn_nodes <- c(pns_num_conn_nodes, sum(rowSums(pns_agg_mx > prop_base_thresh)))
  
  clean_agg_mx <- readRDS(paste0(stblfolder_rd1, 'clean_lambda', lambda, '_numsmps', numsmps, '_propsubsmp0.8.rds'))
  clean_agg_prop_mx <- clean_agg_mx/numsmps
  clean_num_stable_edges <- c(clean_num_stable_edges, sum(clean_agg_prop_mx[upper.tri(clean_agg_prop_mx)] > prop_thresh))
  clean_num_overall_edges <- c(clean_num_overall_edges, sum(clean_agg_prop_mx[upper.tri(clean_agg_prop_mx)] > prop_base_thresh))
  clean_num_conn_nodes <- c(clean_num_stable_edges, sum(rowSums(clean_agg_mx > prop_base_thresh)))
}
metrics_rd1$pns_propStblEdges <- pns_num_stable_edges/pns_num_overall_edges
metrics_rd1$clean_propStblEdges <- clean_num_stable_edges/clean_num_overall_edges

saveRDS(metrics_rd1, paste0(folder_rd1, 'full_metrics_rd1.rds'))

metrics_rd1 <- readRDS(paste0(folder_rd1, 'full_metrics_rd1.rds'))
pdf(paste0(figpath, 'metrics_rd1_smoothed.pdf'), height=10, width=10)
par(mfrow = c(3, 3), oma=c(0,0,3,0))
create_legend <- T
for (i in 1:length(metrics)) {
  metric <- metrics[i]
  metric_title <- metric_titles[i]
  
  if (paste0('pns_', metric) %in% names(metrics_rd1)) {
    pns_vals <- metrics_rd1[[paste0('pns_', metric)]]
    clean_vals <- metrics_rd1[[paste0('clean_', metric)]]
    
    ymin <- min(c(pns_vals, clean_vals), na.rm=T)
    ymax <- max(c(pns_vals, clean_vals), na.rm=T)
    
    plot(lambdas_rd1, clean_vals, xlab='Round 1 Lambdas', ylab=metric_title, main=metric_title, col='blue', ylim=c(ymin, ymax))
    grid(nx = NA, ny = NULL,
         lty = 3,      # Grid line typeb. n 
         col = "lightgray", # Grid line color
         lwd = 1)      # Grid line width
    
    abline(v=lambdas_rd1, col='lightgray', lty=3, lwd=1)
    lines(lambdas_rd1, clean_vals, col='blue')  
    
    
    points(lambdas_rd1, pns_vals, ylab=metric, main=metric, col='black')
    lines(lambdas_rd1, pns_vals, col='black')
    
    print(clean_vals)
    smoothed <- lowess(lambdas_rd1[!is.na(clean_vals)], clean_vals[!is.na(clean_vals)])
    print(smoothed)
    lines(smoothed$x, smoothed$y, col='lightblue')
    
    abline(v=0.15, col='darkred', lty=2)
    abline(v=0.3, col='darkred', lty=2)
    
    if (create_legend) {
      legend("topright", 
             legend = c("PNS", "Clean", "Clean LOWESS"), 
             col = c('black', 'blue', 'lightblue'),
             lty = c(1, 1, 1))
      create_legend <- F
    }
  }
}
mtext("ASD: Round 1 Metrics", side = 3, line = 0.5, outer = TRUE)
dev.off()

## Round 2 ----
folder_rd2 <- 'pns_metrics_rd2/'
dir.create(folder_rd2)
vizfolder_rd2 <- paste0(folder_rd2, 'viz/')
dir.create(vizfolder_rd2)
stblfolder_rd2 <- paste0(folder_rd2, 'stability/')
dir.create(stblfolder_rd2)

lambdas_rd2 <-  seq(0.15, 0.3, 0.01)

### Easier Metrics ----
metrics_rd2 <- construct_graphs_metrics(X, pvals, pvalthresh=NULL, corthresh=0.5, 
                                        pvalscreen=highsubselgenes, lambdas=lambdas_rd2, 
                                        degcut=1, maxnumcomp=NULL, mincompsize=10, degdivide=5,
                                        standardize=T, savefile=folder_rd2, saveviz=vizfolder_rd2, verbose=T)
saveRDS(metrics_rd2, paste0(folder_rd2, 'metrics_rd2.rds'))
metrics_rd2 <- readRDS(paste0(folder_rd2, 'metrics_rd2.rds'))

metrics_rd2$pns_normMeanShortestPath <- metrics_rd2[[paste0('pns_', 'meanShortestPath')]]/metrics_rd2[[paste0('pns_', 'numNodes')]]
metrics_rd2$clean_normMeanShortestPath <- metrics_rd2[[paste0('clean_', 'meanShortestPath')]]/metrics_rd2[[paste0('clean_', 'numNodes')]]

### CV Mean Squared Error ----
# Resulting matrix: (# core genes) x (# lambdas)
# PNS: For lambda, for each core gene, run linear regression on other core genes + nbr genes --> error
# MSE_CV: For each lambda, for each core gene, PNS with K-fold CV --> K-Fold CV error for each gene, average over K folds 
# --> Return matrix of mean K-fold CVs, one for each core gene and one for each lambda
t0 <- Sys.time()
cvmse_mx_rd2 <- get_MSEs_CV(X, pvals, pvalthresh=NULL, corthresh=0.5, pvalscreen=highsubselgenes, lambdas=lambdas_rd2, 
                            standardize=T, verbose=T, K=5)
t1 <- Sys.time()
print(t1-t0)
saveRDS(cvmse_mx_rd2, paste0(folder_rd2, 'cvmse_mx_rd2.rds'))

# Save mean and median CV MSE as both pns and clean for easier plotting
mean_cvmse_rd2 <- colMeans(cvmse_mx_rd2)
metrics_rd2$pns_meanCVMSE <- mean_cvmse_rd2
metrics_rd2$clean_meanCVMSE <- mean_cvmse_rd2

med_cvmse_rd2 <- apply(cvmse_mx_rd2, 2, median)
metrics_rd2$pns_medCVMSE <- med_cvmse_rd2
metrics_rd2$clean_medCVMSE <- med_cvmse_rd2

### Stability ----
# Can run on server
numsmps <- 100
prop_thresh <- 0.6
prop_base_thresh <- 0.1
pns_num_stable_edges <- c()
pns_num_conn_nodes <- c()
pns_num_overall_edges <- c()
clean_num_stable_edges <- c()
clean_num_conn_nodes <- c()
clean_num_overall_edges <- c()
for (i in 1:length(lambdas_rd2)) {
  lambda <- lambdas_rd2[i]
  cat(lambda, '')
  pns_agg_mx <- readRDS(paste0(stblfolder_rd2, 'pns_lambda', lambda, '_numsmps', numsmps, '_propsubsmp0.8.rds'))
  pns_agg_prop_mx <- pns_agg_mx/numsmps
  pns_num_stable_edges <- c(pns_num_stable_edges, sum(pns_agg_prop_mx[upper.tri(pns_agg_prop_mx)] > prop_thresh))
  pns_num_overall_edges <- c(pns_num_overall_edges, sum(pns_agg_prop_mx[upper.tri(pns_agg_prop_mx)] > prop_base_thresh))
  pns_num_conn_nodes <- c(pns_num_conn_nodes, sum(rowSums(pns_agg_mx > prop_base_thresh)))
  
  clean_agg_mx <- readRDS(paste0(stblfolder_rd2, 'clean_lambda', lambda, '_numsmps', numsmps, '_propsubsmp0.8.rds'))
  clean_agg_prop_mx <- clean_agg_mx/numsmps
  clean_num_stable_edges <- c(clean_num_stable_edges, sum(clean_agg_prop_mx[upper.tri(clean_agg_prop_mx)] > prop_thresh))
  clean_num_overall_edges <- c(clean_num_overall_edges, sum(clean_agg_prop_mx[upper.tri(clean_agg_prop_mx)] > prop_base_thresh))
  clean_num_conn_nodes <- c(clean_num_stable_edges, sum(rowSums(clean_agg_mx > prop_base_thresh)))
}
metrics_rd2$pns_propStblEdges <- pns_num_stable_edges/pns_num_overall_edges
metrics_rd2$clean_propStblEdges <- clean_num_stable_edges/clean_num_overall_edges

saveRDS(metrics_rd2, paste0(folder_rd2, 'full_metrics_rd2.rds'))

metrics_rd2 <- readRDS(paste0(folder_rd2, 'full_metrics_rd2.rds'))
pdf(paste0(figpath, 'metrics_rd2_smoothed.pdf'), height=10, width=10)
par(mfrow = c(3, 3), oma=c(0,0,3,0))
create_legend <- T
for (i in 1:length(metrics)) {
  metric <- metrics[i]
  metric_title <- metric_titles[i]
  
  if (paste0('pns_', metric) %in% names(metrics_rd2)) {
    pns_vals <- metrics_rd2[[paste0('pns_', metric)]]
    clean_vals <- metrics_rd2[[paste0('clean_', metric)]]
    
    ymin <- min(c(pns_vals, clean_vals), na.rm=T)
    ymax <- max(c(pns_vals, clean_vals), na.rm=T)
    
    plot(lambdas_rd2, clean_vals, xlab='Round 2 Lambdas', ylab=metric_title, main=metric_title, col='blue', ylim=c(ymin, ymax))
    grid(nx = NA, ny = NULL,
         lty = 3,      # Grid line type
         col = "lightgray", # Grid line color
         lwd = 1)      # Grid line width
    
    abline(v=lambdas_rd2, col='lightgray', lty=3, lwd=1)
    lines(lambdas_rd2, clean_vals, col='blue')  
    
    
    points(lambdas_rd2, pns_vals, ylab=metric, main=metric, col='black')
    lines(lambdas_rd2, pns_vals, col='black')
  
    
    lines(lambdas_rd2, lowess(lambdas_rd2, clean_vals)$y, col='lightblue')
    
    abline(v=0.22, col='red', lty=2)

    if (create_legend) {
      legend("topright", 
             legend = c("PNS", "Clean", "Clean LOWESS"), 
             col = c('black', 'blue', 'lightblue'),
             lty = c(1, 1, 1))
      create_legend <- F
    }
  }
}
mtext("ASD: Round 2 Metrics", side = 3, line = 0.5, outer = TRUE)
dev.off()

for (i in 1:length(lambdas_rd1)) {
  lambda <- lambdas_rd1[i]
  cleangraph <- readRDS(paste0(folder_rd1, 'clean_lambda', lambda, '.rds'))
  saveviz_lambda <- paste0(folder_rd1, 'viz/', i,  '_clean_lambda', lambda, '.pdf')
  set.seed(42)
  if (length(V(cleangraph$net)) > 1750) {
    layout <- layout_with_drl(cleangraph$net)
  } else {
    layout <- layout_with_fr(cleangraph$net)
  }
  pdf(saveviz_lambda, width=10, height=10)
  par(mar=c(0, 0, 0, 0) + .1, 
      bg = "white",  # Set the entire canvas background to white
      fg = "white",  # Set the foreground/border color to white
      bty = "n")     # Suppress the plot box border (if any)
  plot.igraph(cleangraph$net,
              layout = layout,
              vertex.label=NA)
  if (!is.null(dev.list())) {
    dev.off()
  }
}

for (i in 1:length(lambdas_rd2)) {
  lambda <- lambdas_rd2[i]
  cleangraph <- readRDS(paste0(folder_rd2, 'clean_lambda', lambda, '.rds'))
  saveviz_lambda <- paste0(folder_rd2, 'viz/', i,  '_clean_lambda', lambda, '.pdf')
  set.seed(42)
  if (length(V(cleangraph$net)) > 1750) {
    layout <- layout_with_drl(cleangraph$net)
  } else {
    layout <- layout_with_fr(cleangraph$net)
  }
  pdf(saveviz_lambda, width=10, height=10)
  par(mar=c(0, 0, 0, 0) + .1, 
      bg = "white",  # Set the entire canvas background to white
      fg = "white",  # Set the foreground/border color to white
      bty = "n")     # Suppress the plot box border (if any)
  plot.igraph(cleangraph$net,
              layout = layout,
              vertex.label=NA)
  if (!is.null(dev.list())) {
    dev.off()
  }
}

# Select lambda based on metrics
best.lambda <- 0.22

dir.create('final/')
best.params <- list()
best.params$lambda <- best.lambda
saveRDS(best.params, 'final/bestparams.rds')

# Construct and save final graph
pnsgraph <- construct_graph(X, pvals, pvalthresh = NULL, corthresh = 0.5, pvalscreen = highsubselgenes, 
                            lambda = best.lambda, savefile=paste0('final/pns_lambda',best.lambda, '.rds'), standardize=T, verbose=T)
cleangraph <- clean_graph(pnsgraph$graphfinal, degcut=1, mincompsize=10, degdivide=5, savefile=paste0('final/clean_lambda',best.lambda, '.rds'), verbose=T)

set.seed(42)
layout <- layout_with_fr(cleangraph$net)
png(paste0(figpath, 'clean_lambda', best.lambda, '.png'), width=2400, height=2400, res=300)
par(mar=c(0, 0, 0, 0)+.1)
plot.igraph(cleangraph$net,
            layout = layout,
            vertex.label=NA)
dev.off()
