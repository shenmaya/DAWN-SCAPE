# Load Libraries ----
library(here)

# Set Working Directory ----
setwd(here("ADpipeline"))

sapply(paste0(here(), '/mycode/', list.files(paste0(here(), '/mycode/'))), source)

# Set Figure Save Path
figpath <- 'figs/3_RunJointHMRF/'
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

# Load z-scores 
gen_df <- read.csv('data/gen_pvals.csv', row.names=1)
de_df <- read.csv('data/de_pvals.csv', row.names=1)
rownames(gen_df) <- gen_df$ensembl
rownames(de_df) <- de_df$ensembl

head(gen_df)
head(de_df)

#'#######################################################################
# Select Z-Scores ----
#'#######################################################################

gen_zscores <- gen_df[genes,]$zscores
de_zscores <- de_df[genes,]$zscores_cal
names(gen_zscores) <- genes
names(de_zscores) <- genes

plot(gen_zscores, de_zscores)  
dev.off()

#'#######################################################################
# Set Joint-HMRF Hyperparameters  ----
#'#######################################################################

## Define Threshold Hyperparameter Grid ----
pval_seq <- c(0.01, seq(0.05, 0.5, 0.05))
initthres_deseq <- qnorm(1-pval_seq)
initthres_genseq <- qnorm(1-pval_seq)
trimthres_seq <- c(4) 

seq_df <- expand.grid(initthres_deseq, initthres_genseq, trimthres_seq, trimthres_seq)
colnames(seq_df) <- c('initthres_de', 'initthres_gen', 'trimthres_de', 'trimthres_gen')
seq_df

viz_startmx_pval(de_zscores, gen_zscores, pval_seq, pval_seq, 'AD',
                 paste0(figpath, 'lambda', best.lambda, 'start_numactive'), height=2000)

param_numiter <- 100
param_contrun <- F
folder_grid <- 'grid/'
dir.create(folder_grid)
folder_gridfiles <- paste0(folder_grid, 'load/')
dir.create(folder_gridfiles)

## Create list of parameters ----
# Cycle through 2D grid and check if any parameter set violate setup
params_list <- list()
ctr <- 1
num_active <- c()
init_state1_genes <- list()
init_state2_genes <- list()
init_state3_genes <- list()
init_state4_genes <- list()
for (i in 1:nrow(seq_df)) {
  
  params <- seq_df[i, ]
  
  initthres_de <- params$initthres_de
  initthres_gen <- params$initthres_gen
  trimthres_de <- params$trimthres_de
  trimthres_gen <- params$trimthres_gen
  
  Istart <- 1*((de_zscores > initthres_de) & (gen_zscores > initthres_gen)) + 2*((de_zscores > initthres_de) & (gen_zscores <= initthres_gen)) + 3*((de_zscores <= initthres_de) & (gen_zscores > initthres_gen)) + 4*((de_zscores <= initthres_de) & (gen_zscores <= initthres_gen))
  fixindex <- 1*((de_zscores > trimthres_de) & (gen_zscores > trimthres_gen)) + 2*((de_zscores > trimthres_de) & (gen_zscores <= trimthres_gen))+ 3*(de_zscores <= trimthres_de) 
  
  num_active <- c(num_active, sum(Istart == 1))
  
  init_state1_genes[[i]] <- names(which(Istart == 1))
  init_state2_genes[[i]] <- names(which(Istart == 2))
  init_state3_genes[[i]] <- names(which(Istart == 3))
  init_state4_genes[[i]] <- names(which(Istart == 4))
  
  mu1 <- mean(de_zscores[Istart <= 2 & fixindex == 3]) 
  sigma1 <- (sd(de_zscores[Istart <= 2]))^2
  
  mu2 <- mean(gen_zscores[(Istart %% 2==1) & fixindex == 3])
  sigma2 <- (sd(gen_zscores[(Istart %% 2==1)]))^2
  
  sigma0 <- (sd(c(gen_zscores[(Istart %% 2)==0], de_zscores[Istart > 2])))^2
  
  if (any(is.na(c(mu1, sigma1, mu2, sigma2, sigma0)))) {
    print('starting NA')
    print(c(mu1, sigma1, mu2, sigma2, sigma0))
  } else {
    params_list[[ctr]] <- list(initthres_de, initthres_gen, # initthres
                               trimthres_de, trimthres_gen, # trimthres
                               param_numiter, # numiter
                               param_contrun, # contrun
                               paste0(folder_grid, 'savejointhmrf_', i, '.rds'),
                               paste0(folder_grid, 'graphviz_', i, '.pdf'))
    ctr <- ctr + 1
  }
}
init_state_genes <- list()
init_state_genes[[1]] <- init_state1_genes
init_state_genes[[2]] <- init_state2_genes
init_state_genes[[3]] <- init_state3_genes
init_state_genes[[4]] <- init_state4_genes
saveRDS(init_state_genes, paste0(folder_grid, 'init_state_genes.rds'))

# Save list of parameters
saveRDS(params_list, paste0(folder_gridfiles, 'params_list_', param_numiter, 'iters.rds'))
saveRDS(seq_df, paste0(folder_gridfiles, 'seq_df_', param_numiter, 'iters.rds'))

zscores <- list()
zscores$de <- de_zscores
zscores$gen <- gen_zscores
saveRDS(zscores, paste0(folder_gridfiles, 'zscores.rds'))

set.seed(42)
# layout <- layout_with_fr(cleanraw$net)
layout <- layout_with_drl(cleanraw$net)
saveRDS(cleanraw$net, paste0(folder_gridfiles, 'graph_lambda', best.lambda, '.rds'))
saveRDS(layout, paste0(folder_gridfiles, 'layout_seed42.rds'))

print(sort(zscores$de, decreasing=T)[1:20])
print(sort(zscores$gen, decreasing=T)[1:20])
cor(zscores$de, zscores$gen)

#'#######################################################################
# Run Joint-HMRF ----
#'#######################################################################

# Recommend running next part on server
library(igraph)
library(fdrtool)
library(splines)
library(parallel)
library(mixR)
library(here)

setwd(here("ADpipeline"))

sapply(paste0(here(), '/mycode/', list.files(paste0(here(), '/mycode/'))), source)

folder <- paste0('grid/')
numiters <- 100
folder_load <- paste0(folder, 'load/')

params_list <- readRDS(paste0(folder_load, 'params_list_', numiters, 'iters.rds'))
seq_df <- readRDS(paste0(folder_load, 'seq_df_', numiters, 'iters.rds'))

zscores <- readRDS(paste0(folder_load, 'zscores.rds'))
zscorede <- zscores$de
zscoregen <- zscores$gen
print(sort(zscorede, decreasing=T)[1:20])
print(sort(zscoregen, decreasing=T)[1:20])
cor(zscores$de, zscores$gen)

graph <- readRDS(paste0(folder_load, 'graph_lambda0.16.rds'))
layout <- readRDS(paste0(folder_load, 'layout_seed42.rds'))

## Run Parallelized Joint-HMRF ----
detectCores()
# Set number of cores 
num_cores <- 50

# Toy test
dir.create(paste0(folder, 'toy/'))
test <- params_list[[100]]
test[[5]] <- 2 # 2 iters
test[[7]] <- paste0(folder, 'toy/savejointhmrf_100_test.rds') # Modify save str
test[[8]] <- paste0(folder, 'toy/savejointhmrf_100_test.pdf') # Modify save fig str
run_hmrf_parallel(test)

t0 <- Sys.time()
result <- mclapply(params_list, run_hmrf_parallel, mc.cores = num_cores)
t1 <- Sys.time()
print(t1-t0)

#'#######################################################################
# Analyze Joint-HMRF Results ----
#'#######################################################################
# Back to local
folder_grid <- 'grid/'
folder_gridfiles <- paste0(folder_grid, 'load/')
param_numiter <- 100

params_list <- readRDS(paste0(folder_gridfiles, 'params_list_', param_numiter, 'iters.rds'))
seq_df <- readRDS(paste0(folder_gridfiles, 'seq_df_', param_numiter, 'iters.rds'))
init_state_genes <- readRDS(paste0(folder_grid, 'init_state_genes.rds'))
init_state1_genes <- init_state_genes[[1]]
init_state2_genes <- init_state_genes[[2]] 
init_state3_genes <- init_state_genes[[3]] 
init_state4_genes <- init_state_genes[[4]] 

savefigpath <- paste0(figpath, folder_grid)
dir.create(savefigpath)

seq_df$p_de <- 1-pnorm(seq_df$initthres_de)
seq_df$p_gen <- 1-pnorm(seq_df$initthres_gen)

pval_seq <- round(sort(unique(1-pnorm(seq_df$initthres_de))), 2)

num_active_genes <- c()
num_reactive_genes <- c()
num_state3_genes <- c()
finallogpposts <- c()
param_bound <- c()
active_genes_list <- list()
states_list <- list()
b01s <- c()
b02s <- c()
b03s <- c()
b11s <- c()
b12s <- c()
mu1s <- c()
mu2s <- c()
sigma1s <- c()
sigma2s <- c()
sigma01s <- c()
sigma02s <- c()
num_active_overlap <- c()
for (i in 1:length(params_list)) {
  hmrf_res <- readRDS(paste0(params_list[[i]][[7]]))
  
  # There is some oscillation, take iteration t = 100 vs 99 based on which has larger pseudo-posterior
  if (hmrf_res$logppost_hist[100] >= hmrf_res$logppost_hist[99]) {
    t <- 101
  } else {
    t <- 100
  }
  
  hrmf_final <- final_states(genes, as_adjacency_matrix(cleanraw$net), hmrf_res$state_hist[t,], 
                             de_zscores, gen_zscores,
                             hmrf_res$param_hist[t, 1], hmrf_res$param_hist[t, 2], hmrf_res$param_hist[t, 3], 
                             hmrf_res$param_hist[t, 4], hmrf_res$param_hist[t, 5], 
                             hmrf_res$param_hist[t, 6], hmrf_res$param_hist[t, 7], hmrf_res$param_hist[t, 8], hmrf_res$param_hist[t, 9],
                             hmrf_res$param_hist[t, 10], hmrf_res$param_hist[t, 11])$Iupdate
  names(hrmf_final) <- genes
  hmrf_res$Iupdate <- hrmf_final
  
  states_list[[i]] <- hmrf_res$Iupdate
  
  num_active <- sum(hmrf_res$Iupdate == 1)
  num_active_genes <- c(num_active_genes, num_active)
  active_genes_list[[i]] <- names(hmrf_res$Iupdate)[which(hmrf_res$Iupdate == 1)]
  
  num_reactive <- sum(hmrf_res$Iupdate == 2)
  num_reactive_genes <- c(num_reactive_genes, num_reactive)
  
  num_state3_genes <- c(num_state3_genes, sum(hmrf_res$Iupdate == 3))
  num_active_overlap <- c(num_active_overlap, length(intersect(init_state1_genes[[i]], names(which(hmrf_res$Iupdate == 1)))))
  
  logppost <- hmrf_res$logppost_hist[t-1]
  finallogpposts <- c(finallogpposts, logppost)
  
  b01s <- c(b01s, hmrf_res$param_hist[t, 1])
  b02s <- c(b02s, hmrf_res$param_hist[t, 2])
  b03s <- c(b03s, hmrf_res$param_hist[t, 3])
  b11s <- c(b11s, hmrf_res$param_hist[t, 4])
  b12s <- c(b12s, hmrf_res$param_hist[t, 5])
  mu1s <- c(mu1s, hmrf_res$param_hist[t, 6])
  mu2s <- c(mu2s, hmrf_res$param_hist[t, 8])
  sigma1s <- c(sigma1s, hmrf_res$param_hist[t, 7])
  sigma2s <- c(sigma2s, hmrf_res$param_hist[t, 9])
  sigma01s <- c(sigma01s, hmrf_res$param_hist[t, 10])
  sigma02s <- c(sigma02s, hmrf_res$param_hist[t, 11])
 
  # Check for parameter values near bounds which indicates poor fit, take average over last 5 iterations (from t)
  # For the AD application, we believe very few or no state 3 values is possible and, thus, accept an extreme b03
  if (any(colMeans(hmrf_res$param_hist[(t-5):t, 1:2]) < (-20 + 1e-6)) | any(colMeans(hmrf_res$param_hist[(t-5):t, 4:5]) > (10 - 1e-6))) {
    param_bound <- c(param_bound, NA)
  } else {
    param_bound <- c(param_bound, 1)
  }
  
  print(paste0(i, ' (iter ', t-1, ')'))
  print(tail(hmrf_res$logppost_hist, 2))
  print(hmrf_res$param_hist[(t-5):t,])
  print(table(hmrf_res$Iupdate))
  print('==================')
}

# Which runs should perhaps be ignored due to hidden parameters that are too close to bounds
na_mx <- matrix(param_bound, nrow = length(pval_seq), ncol = length(pval_seq), byrow = TRUE)
na_mx

## Number of active and p-reactive states ----
numactive_mx <- matrix(num_active_genes, nrow = length(pval_seq), ncol = length(pval_seq), byrow = TRUE)
rownames(numactive_mx) <- paste0('gen_', pval_seq)
colnames(numactive_mx) <- paste0('de_', pval_seq)
viz_mx_pval(val_mx=numactive_mx, 
            pval_seq_de=pval_seq, pval_seq_gen=pval_seq,
            title=paste0('AD: Num Active'), saveviz=paste0(savefigpath, "numactive.jpg"))
viz_mx_pval(val_mx=numactive_mx, 
            pval_seq_de=pval_seq, pval_seq_gen=pval_seq, null_mx = na_mx,
            title=paste0('AD: Num Active'), saveviz=paste0(savefigpath, "numactive_na.jpg"))

numreactive_mx <- matrix(num_reactive_genes, nrow = length(pval_seq), ncol = length(pval_seq), byrow = TRUE)
rownames(numreactive_mx) <- paste0('gen_', pval_seq)
colnames(numreactive_mx) <- paste0('de_', pval_seq)
viz_mx_pval(val_mx=numreactive_mx, 
            pval_seq_de=pval_seq, pval_seq_gen=pval_seq,
            title=paste0('AD: Num Reactive'), saveviz=paste0(savefigpath, "numreactive.jpg"))
viz_mx_pval(val_mx=numreactive_mx, 
            pval_seq_de=pval_seq, pval_seq_gen=pval_seq, null_mx = na_mx,
            title=paste0('AD: Num Reactive'), saveviz=paste0(savefigpath, "numreactive_na.jpg"))

## Parameter convergence ----
viz_paramconv_iter(params_list, 'AD: 100 iters', paste0(savefigpath, "paramconv_100iters"), numiters=NULL)

## Active Gene Stability ----
# Average number of stable active genes
n_grid <- length(pval_seq)
index_matrix <- matrix(1:(length(pval_seq)**2), nrow = n_grid, ncol = n_grid, byrow = TRUE)

# Prepare a matrix to store the average overlap with neighbors
overlap_mx <- matrix(NA, nrow = n_grid, ncol = n_grid)

# Loop over grid
for (row in 1:n_grid) {
  for (col in 1:n_grid) {
    
    this_index <- index_matrix[row, col]
    this_genes <- active_genes_list[[this_index]]
    
    neighbor_indices <- c()
    if (row > 1) neighbor_indices <- c(neighbor_indices, index_matrix[row - 1, col])  # Up
    if (row < n_grid) neighbor_indices <- c(neighbor_indices, index_matrix[row + 1, col])  # Down
    if (col > 1) neighbor_indices <- c(neighbor_indices, index_matrix[row, col - 1])  # Left
    if (col < n_grid) neighbor_indices <- c(neighbor_indices, index_matrix[row, col + 1])  # Right
    
    overlaps <- sapply(neighbor_indices, function(neigh_idx) {
      length(intersect(this_genes, active_genes_list[[neigh_idx]]))
    })
    
    # Option: mean overlap with neighbors
    overlap_mx[row, col] <- mean(overlaps)
    
    # If you want sum instead, use:
    # overlap_mx[row, col] <- sum(overlaps)
  }
}
rownames(overlap_mx) <- paste0('gen_', pval_seq)
colnames(overlap_mx) <- paste0('de_', pval_seq)
viz_mx_pval(val_mx=overlap_mx, 
            pval_seq_de=pval_seq, pval_seq_gen=pval_seq,
            title=paste0('AD: Num Active Overlapping'), saveviz=paste0(savefigpath, "overlap_numactive.jpg"))
viz_mx_pval(val_mx=overlap_mx, 
            pval_seq_de=pval_seq, pval_seq_gen=pval_seq,
            title=paste0('AD: Num Active Overlapping'), saveviz=paste0(savefigpath, "overlap_numactive_na.jpg"),
            null_mx = na_mx)

# Proportion of stable active genes
# Prepare a matrix to store the average overlap with neighbors
overlapprop_mx <- matrix(NA, nrow = n_grid, ncol = n_grid)

# Loop over grid
for (row in 1:n_grid) {
  for (col in 1:n_grid) {
    
    this_index <- index_matrix[row, col]
    this_genes <- active_genes_list[[this_index]]
    
    neighbor_indices <- c()
    if (row > 1) neighbor_indices <- c(neighbor_indices, index_matrix[row - 1, col])  # Up
    if (row < n_grid) neighbor_indices <- c(neighbor_indices, index_matrix[row + 1, col])  # Down
    if (col > 1) neighbor_indices <- c(neighbor_indices, index_matrix[row, col - 1])  # Left
    if (col < n_grid) neighbor_indices <- c(neighbor_indices, index_matrix[row, col + 1])  # Right
    
    nbr_genes <- unique(unlist(active_genes_list[neighbor_indices]))
    
    # Option: mean overlap with neighbors
    overlapprop_mx[row, col] <- length(intersect(this_genes, nbr_genes))/length(nbr_genes)
    
    # If you want sum instead, use:
    # overlap_matrix[row, col] <- sum(overlaps)
  }
}
rownames(overlapprop_mx) <- paste0('gen_', pval_seq)
colnames(overlapprop_mx) <- paste0('de_', pval_seq)

viz_mx_pval(val_mx=overlapprop_mx, 
            pval_seq_de=pval_seq, pval_seq_gen=pval_seq,
            title=paste0('AD: Prop Active Overlapping'), saveviz=paste0(savefigpath, "overlap_propactive.jpg"))
viz_mx_pval(val_mx=overlapprop_mx, 
            pval_seq_de=pval_seq, pval_seq_gen=pval_seq,
            title=paste0('AD: Prop Active Overlapping'), saveviz=paste0(savefigpath, "overlap_propactive_na.jpg"),
            null_mx = na_mx)

## All state stability: Hamming distance ----
# Prepare a matrix to store average Hamming distance with neighbors
hamdist_mx <- matrix(NA, nrow = n_grid, ncol = n_grid)

# Loop over grid
for (row in 1:n_grid) {
  for (col in 1:n_grid) {
    
    this_index <- index_matrix[row, col]
    this_states <- states_list[[this_index]]
    
    neighbor_indices <- c()
    if (row > 1) neighbor_indices <- c(neighbor_indices, index_matrix[row - 1, col])  # Up
    if (row < n_grid) neighbor_indices <- c(neighbor_indices, index_matrix[row + 1, col])  # Down
    if (col > 1) neighbor_indices <- c(neighbor_indices, index_matrix[row, col - 1])  # Left
    if (col < n_grid) neighbor_indices <- c(neighbor_indices, index_matrix[row, col + 1])  # Right
    
    nbr_hamdists <- c()
    for (nbr_index in neighbor_indices) {
      nbr_states <- states_list[[nbr_index]]
      nbr_hamdists <- c(nbr_hamdists, abs(this_states - nbr_states))
    }
    
    hamdist_mx[row, col] <- mean(nbr_hamdists)
  }
}
rownames(hamdist_mx) <- paste0('gen_', pval_seq)
colnames(hamdist_mx) <- paste0('de_', pval_seq)

viz_mx_pval(val_mx=hamdist_mx, 
            pval_seq_de=pval_seq, pval_seq_gen=pval_seq,
            title=paste0('AD: Avg Hamm Dist'), saveviz=paste0(savefigpath, "overlap_hamdist.jpg"))
viz_mx_pval(val_mx=hamdist_mx, 
            pval_seq_de=pval_seq, pval_seq_gen=pval_seq,
            title=paste0('AD: Avg Hamm Dist'), saveviz=paste0(savefigpath, "overlap_hamdist_na.jpg"),
            null_mx = na_mx)

# Don't consider outer layer due to boundary effect 
rownames(index_matrix) = colnames(index_matrix) = pval_seq
index_matrix

# Select Final Runs ----
# Selected final run idxs based on Num Active, Prop Overlapping (Active), and Hamming Dist (All State Types) 
dir.create('final/states/')
final_idxs <- c(96, 97, 107, 108, 109)
saveRDS(final_idxs, paste0('final/states/final_idxs.rds'))

# Save final states from these runs
final_activegenes <- list()
final_genes <- list()
for (i in 1:length(final_idxs)) {
  idx <- final_idxs[i]
  states <- states_list[[idx]]
  final_activegenes[[paste0('idx', idx)]] <- names(which(states == 1))
  final_genes[[paste0('idx', idx)]] <- states
}
lapply(final_activegenes, length)
length(Reduce(union, final_activegenes))
length(Reduce(intersect, final_activegenes))

saveRDS(final_activegenes, paste0('final/states/final_activegenes.rds'))
saveRDS(final_genes, paste0('final/states/final_genes.rds'))
