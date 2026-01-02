# Functions that should be in this file:
# graph construction main function
# pns algorithm
# clean graph function 

library(glmnet)
library(igraph)
library(doParallel)

#' Construct network using p-value and correlation thresholding and PNS (Alg 1)
#' 
#' @param X Gene expression matrix, samples x genes
#' @param pvals Vector of p-values
#' @param pvalthresh Gene p-value threshold
#' @param corthresh Neighbor correlation threshold
#' @param pvalscreen Vector of gene names corresponding to p-value thresholding (screen1_genes), default=NULL
#' @param lambda Scalar (CV=F) or vector (CV=T) of regularization parameters, 
#' @param standardize Whether or not to standardize data for PNS, default=T
#' @param savefile String denoting full path and file name for where to store file (rds), default=NULL
#' @param verbose Whether or not to print out progress messages, default=T
#' @returns If CV=F, a list containing
#' - genesfinal: Vector of final genes in network
#' - graphfinal: 0-1 Matrix denoting final gene network
#' - coef_graphfinal: Weighted matrix with coef values as edges 
#' If CV=T, vector of MSEs from K-fold CV across lambdas
construct_graph <- function(X, pvals, pvalthresh, corthresh, pvalscreen, lambda, 
                            standardize=T, savefile=NULL, verbose=T) {
  
  if (!is.null(savefile)) {
    savefilesplit <- strsplit(savefile, '\\.')[[1]]
    if (tail(savefilesplit, 1) != 'rds') {
      warning('`savefile` does not end in .rds, manually appending .rds.')
      savefile <- paste0(savefile, '.rds')
    }
  }
  
  # P-val screening and neighbor retrieval (Alg 1, Steps 1-3)
  pns_inputs <- screen_retrieve(X, pvals, pvalthresh, corthresh, pvalscreen, verbose=T)
  
  # Construct graph via PNS (Alg 1, Step 4)
  final <- run_pns(pns_inputs$X, pns_inputs$core_genes, pns_inputs$nbr_genes, lambda, 
                     standardize=standardize, savefile=savefile, verbose=verbose)
  
  return(final)
}

#' Conduct p-value and correlation screening and retrieve neighbors (Alg 1, Steps 1-3)
#' 
#' @param X Gene expression matrix, samples x genes
#' @param pvals Vector of p-values
#' @param pvalthresh Gene p-value threshold
#' @param corthresh Neighbor correlation threshold
#' @param pvalscreen Vector of gene names corresponding to p-value thresholding (screen1_genes), default=NULL
#' @param standardize Whether or not to standardize data for PNS, default=T
#' @param verbose Whether or not to print out progress messages, default=T
#' @returns A list containing inputs to run_pns:
#' - X: Gene expression matrix of genes in set V, samples x genes
#' - core_genes: Vector containing the genes in set S
#' - nbr_genes: Vector containing the genes in set V / S
screen_retrieve <- function(X, pvals, pvalthresh, corthresh, pvalscreen, verbose=T) {
  
  if (is.null(colnames(X)) | is.null(names(pvals))) {
    stop('Please provide colnames for `X` and names for `pvals`')
  }
  if (dim(X)[2] != length(pvals)) {
    stop('Number of genes in `X` and `pvals` do not match')
  }
  if (length(unique(colnames(X), names(pvals))) != length(names(pvals))) {
    stop('Genes in `X` (colnames) and `pvals` do not match')
  }
  if (!is.null(pvalscreen) & !is.null(pvalthresh)) {
    warning('Both `pvalscreen` and `pvalthresh` inputted, using `pvalscreen`')
  }
  if (!is.null(pvalthresh)) {
    if ((pvalthresh > 1) | (pvalthresh < 0)) {
      stop('Invalid `pvalthresh`, must be in [0, 1]')
    }
  }
  
  # Make sure X and pvals have the same gene order
  genes <- names(pvals)
  X <- X[, genes]
  if (verbose) {
    cat('Number of starting genes:', length(genes), '\n')
  }
  
  # P-Value Screening: Exclude genes with p-value < pvalthresh (Alg 1, Step 1)
  if (!is.null(pvalscreen)) {
    screen1_genes <- pvalscreen
  } else {
    screen1_genes <- unique(genes[which(pvals <= pvalthresh)])
  }
  if (verbose) {
    if (!is.null(pvalscreen)) {
      cat('Number of genes not in `pvalscreen`:', length(genes) - length(screen1_genes), '\n')
    } else {
      cat(paste0('Number of genes that do not meet `pvalthresh`=', pvalthresh, ': ', length(genes) - length(screen1_genes)))
    }
    cat(paste0('Number of key genes: ', length(screen1_genes), '\n'))
  }
  
  # Corr Screening (Alg 1, Step 2)
  corX <- abs(cor(X))
  corX <- (corX > corthresh) + 0
  diag(corX) <- 0
  screen1_corX <- corX[screen1_genes, screen1_genes]
  # Exclude isolated nodes
  core_genes <- unique(screen1_genes[which(rowSums(screen1_corX) != 0)])
  if (verbose) {
    cat(paste0('Number of additional genes that do not meet corthresh=', corthresh, ': ', length(screen1_genes) - length(core_genes), '\n'))
    cat(paste0('Number of core genes: ', length(core_genes), '\n'))
  }
  # Neighbor Retrieval (Alg 1, Step 3)
  ret_nbrs <- function(gene_i) {
    return(genes[which(corX[gene_i, ] == 1)])
  }
  nbr_genes <- unique(unlist(lapply(core_genes, ret_nbrs)))
  nbr_genes <- setdiff(nbr_genes, core_genes)
  
  if (verbose) {
    cat(paste0('Number of neighbor genes: ', length(nbr_genes),  '\n'))
    cat(paste0('Number of genes for PNS: ', length(core_genes)+length(nbr_genes), '\n'))
  }
  
  final <- list()
  final$X <- X
  final$core_genes <- core_genes
  final$nbr_genes <- nbr_genes
  return(final)
}

#' Function to run PNS (Alg 1, Step 4)
#' 
#' @param X Gene expression matrix of genes in set V, samples x genes
#' @param core_genes Vector containing the genes in set S
#' @param nbr_genes Vector containing the genes in set V / S
#' @param lambda Regularization parameter
#' @param standardize Whether or not to standardize data, default=T
#' @param savefile String denoting full path and file name for where to store file (rds), default=NULL
#' @param verbose Whether or not to print out progress messages, default=T
#' @returns A list containing
#' - genesfinal: Vector of final genes in network
#' - graphfinal: 0-1 Matrix denoting final gene network
#' - coef_graphfinal: Weighted matrix with coef values as edges 
run_pns <- function(X, core_genes, nbr_genes, lambda, standardize=T, savefile=NULL, verbose=F) {
  
  numcg <- length(core_genes)
  numng <- length(nbr_genes)
  
  if (numcg == 0) {
    stop('No core genes inputted (`core_genes` cannot be empty).')
  }
  if (numng == 0) {
    warning('No neighbor genes inputted (`nbr_genes` is empty).')
  }
  if (length(intersect(core_genes, nbr_genes))) {
    stop('`core_genes` and `nbr_genes` should have no intersection.')
  }
  if (!all(c(core_genes, nbr_genes) %in% colnames(X))) {
    stop('Not all `core_genes` and `nbr_genes` are in `X`.')
  }
  
  genes <- c(core_genes, nbr_genes)
  # Expression data of core genes
  yexp <- X[, core_genes] 
  # Expression data of neighbor genes 
  if (numng > 0) {
    xexp <- X[, nbr_genes] 
    xyexp <- cbind(yexp, xexp)
  } else{
    xyexp <- yexp
  }
  
  edge_coef <- function(i) {
    vectorr <- rep(0, numcg + numng)
    model <- glmnet(xyexp[,-i], 
                     yexp[,i], 
                     family = "gaussian", lambda = lambda, 
                     standardize = standardize)
    vectorr[-i] <- c(model$beta[, 1])
    return(vectorr)
  }
  
  edgecoef_mx <- sapply(1:numcg, edge_coef)
  if (numng > 0) {
    edgecoef_mx <- cbind(edgecoef_mx, matrix(0, nrow=numcg+numng, ncol=numng))
  }
  rownames(edgecoef_mx) = colnames(edgecoef_mx) = genes
  
  # Binarize matrix - E matrix
  binary_mx <- edgecoef_mx
  binary_mx <- (abs(binary_mx) > 0) * 1
  # Symmetrize matrix - Omega matrix
  binary_mx <- 1 - (1 - binary_mx) * (1 - t(binary_mx))
  
  final_genes <- genes #genes[which(rowSums(binary_mx) != 0)]
  final_mx <- binary_mx[final_genes, final_genes]
  finalcoef_mx <- edgecoef_mx[final_genes, final_genes]
  
  if (verbose) {
    cat(paste0('Number of edges among core genes: ', sum(final_mx[core_genes, core_genes])/2, '\n'))
  }
  
  # if (verbose) {
  #   cat(paste0('Number of edges among core genes: ', sum(final_mx[intersect(core_genes, final_genes), intersect(core_genes, final_genes)])/2, '\n'))
  # }
  
  final <- list()
  final$genesfinal <- final_genes
  final$graphfinal <- final_mx
  final$coef_graphfinal <- finalcoef_mx
  
  # Save final object
  if (!is.null(savefile)) {
    saveRDS(final, savefile)
  }
  return(final)
}

#' Get CV Mean Squared Error for PNS fitting (Alg 1)
#' 
#' @param X Gene expression matrix, samples x genes
#' @param pvals Vector of p-values
#' @param pvalthresh Gene p-value threshold
#' @param corthresh Neighbor correlation threshold
#' @param pvalscreen Vector of gene names corresponding to p-value thresholding (screen1_genes), default=NULL
#' @param lambdas Vector of scalars of regularization parameters, 
#' @param standardize Whether or not to standardize data for PNS, default=T
#' @param verbose Whether or not to print out progress messages, default=T
#' @param K What K should be for K-fold CV, ignored if CV=F, default=5
#' @returns Vector of MSEs from K-fold CV across lambdas averaged across genes
get_MSEs_CV <- function(X, pvals, pvalthresh, corthresh, pvalscreen, lambdas, 
                            standardize=T, verbose=T, K=5) {
  
  # P-val screening and neighbor retrieval (Alg 1, Steps 1-3)
  pns_inputs <- screen_retrieve(X, pvals, pvalthresh, corthresh, pvalscreen, verbose=T)
  
  X <- pns_inputs$X
  core_genes <- pns_inputs$core_genes
  nbr_genes <- pns_inputs$nbr_genes
  
  numcg <- length(core_genes)
  numng <- length(nbr_genes)
    
  if (numcg == 0) {
    stop('No core genes inputted (`core_genes` cannot be empty).')
  }
  if (numng == 0) {
    warning('No neighbor genes inputted (`nbr_genes` is empty).')
  }
  if (length(intersect(core_genes, nbr_genes))) {
    stop('`core_genes` and `nbr_genes` should have no intersection.')
  }
  if (!all(c(core_genes, nbr_genes) %in% colnames(X))) {
    stop('Not all `core_genes` and `nbr_genes` are in `X`.')
  }
  
  genes <- c(core_genes, nbr_genes)
  # Expression data of core genes
  yexp <- X[, core_genes] 
  # Expression data of neighbor genes 
  if (numng > 0) {
    xexp <- X[, nbr_genes] 
    xyexp <- cbind(yexp, xexp)
  } else{
    xyexp <- yexp
  }
  
  # Store mean cross-validated errors
  cvms_mx <- matrix(NA, nrow=length(core_genes), ncol=length(lambdas))
  for (i in 1:length(core_genes)) {
    cv_model <- cv.glmnet(x=xyexp[,-i], 
                          y=yexp[,i], 
                          family = "gaussian", lambda=lambdas, 
                          standardize = standardize, nfolds=K)
    cvms_mx[i,] <- cv_model$cvm
  }
  return(cvms_mx)
}

#' Construct network using p-value and correlation thresholding and PNS (Alg 1)
#' 
#' @param X Gene expression matrix, samples x genes
#' @param pvals Vector of p-values
#' @param pvalthresh Gene p-value threshold
#' @param corthresh Neighbor correlation threshold
#' @param pvalscreen Vector of gene names corresponding to p-value thresholding (screen1_genes), default=NULL
#' @param lambdas Vector of scalars of regularization parameters
#' @param numsmps Number of subsamples, default=100
#' @param propsubsmp Proportion of samples to subsample, default=0.8
#' @param num_cores Number of cores to use, default=NULL
#' @param standardize Whether or not to standardize data for PNS, default=T
#' @param savefile String denoting path and file name for where to store file, '_lambda<lambda>.rds' appended at end, default=NULL
#' @param verbose Whether or not to print out progress messages, default=T
#' @returns Vector of MSEs from K-fold CV across lambdas averaged across genes
construct_graphs_stability <- function(X, pvals, pvalthresh, corthresh, pvalscreen, lambdas, numsmps=100, propsubsmp=0.8, 
                                       degcut=1, maxnumcomp=NULL, mincompsize=10, 
                                       num_cores=NULL, savefile=NULL, standardize=T, verbose=T) {
  
  # P-val screening and neighbor retrieval (Alg 1, Steps 1-3)
  pns_inputs <- screen_retrieve(X, pvals, pvalthresh, corthresh, pvalscreen, verbose=verbose)
  
  # for (lambda in lambdas) {
  #   
  #   # Construct graph via PNS (Alg 1, Step 4)
  #   pnsgraph <- run_pns(pns_inputs$X, pns_inputs$core_genes, pns_inputs$nbr_genes, lambda, 
  #                       standardize=standardize, savefile=savefile_lambda, verbose=verbose)
  #   
  #   cleangraph <- clean_graph(pnsgraph$graphfinal, degcut=2, degdivide=5, savefile=paste0(folder,'clean_lambda',lambda, '.rds'), verbose=T)
  #   
  # }
  
  # Parallelize over lambdas, keeping numsmps sequential
  if (is.null(num_cores)) {
    num_cores <- min(parallel::detectCores() - 2, length(lambdas))
  }
  cat('Using', num_cores, 'cores \n')
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  on.exit({
    cat("Stopping cluster...\n")
    parallel::stopCluster(cl)
  }, add = TRUE)
  
  results <- foreach(i = 1:length(lambdas), .packages = c("glmnet", "igraph"), .export = c("run_pns", "clean_graph")) %dopar% {
    cat(i, ' ')
    lambda <- lambdas[i]
    t0 <- Sys.time()
    set.seed(42 + i)  # Ensure different seeds for different lambdas
    
    pns_agg_mx <- NULL  # Initialize pns_agg_mx
    clean_agg_mx <- NULL # Initialize clean_agg_mx
    for (s in 1:numsmps) {  # Run sequentially within each worker
      cat(s)
      smp <- sample(1:nrow(pns_inputs$X), round(nrow(pns_inputs$X) * propsubsmp), replace = FALSE)
      pns_smp_mx <- run_pns(pns_inputs$X[smp, ], pns_inputs$core_genes, pns_inputs$nbr_genes, lambda, 
                        standardize=TRUE, savefile=NULL, verbose=FALSE)
      pns_agg_mx <- if (is.null(pns_agg_mx)) pns_smp_mx$graphfinal else pns_agg_mx + pns_smp_mx$graphfinal
      
      clean_smp_mx <- tryCatch(
        {
          clean_graph(pns_smp_mx$graphfinal, degcut=degcut, maxnumcomp=maxnumcomp, mincompsize=mincompsize, 
                     savefile=NULL, verbose=FALSE)
        }, error = function(msg){
          return(NULL)
        })
      if (is.null(clean_smp_mx)) {
        clean_smp_mx_mod <- matrix(0, nrow=nrow(pns_agg_mx), ncol=ncol(pns_agg_mx))
      } else {
        clean_smp_mx_mod <- matrix(0, nrow=nrow(pns_agg_mx), ncol=ncol(pns_agg_mx))
        rownames(clean_smp_mx_mod) <- colnames(clean_smp_mx_mod) <- rownames(pns_smp_mx$graphfinal)
        clean_smp_mx_mod[rownames(clean_smp_mx$ne), colnames(clean_smp_mx$ne)] <- clean_smp_mx$ne
      }
      clean_agg_mx <- if (is.null(clean_agg_mx)) clean_smp_mx_mod else clean_agg_mx + clean_smp_mx_mod
    }
    t1 <- Sys.time()
    print(t1 - t0)
    
    # Save results
    if (!is.null(savefile)) {
      saveRDS(pns_agg_mx, paste0(savefile, 'pns_lambda', lambda, '_numsmps', numsmps, '_propsubsmp', propsubsmp, '.rds'))
      saveRDS(clean_agg_mx, paste0(savefile, 'clean_lambda', lambda, '_numsmps', numsmps, '_propsubsmp', propsubsmp, '.rds'))
    }
    list(lambda = lambda, pns_agg_mx = pns_agg_mx, clean_agg_mx = clean_agg_mx)
  }
  # stopCluster(cl)
  
  return(results)
}

#' Construct network using p-value and correlation thresholding and PNS (Alg 1)
#' 
#' @param graph Gene expression matrix, samples x genes
#' @param metric String of metric to calculate
#' @returns Metric value of graph 
calc_metric <- function(graph, metric) {
  
  if (!(metric %in% c('R2', 'numEdges', 'numNodes', 'meanMeanEdges', 'meanDegree', 'medianMeanEdges', 'medianDegree', 'meanShortestPath', 'propDeg1'))) {
    stop("Invalid `metric`. Must be one of the following:  c('R2', 'numEdges', 'numNodes', 'meanMeanEdges', 'meanDegree', 'medianMeanEdges', 'medianDegree', 'meanShortestPath', 'propDeg1')")
  }
  
  if (metric == 'R2') {
    return(get_R2(graph))
  } else if (metric == 'numEdges') {
    return(sum(graph)/2)
  } else if (metric == 'numNodes') {
    return(nrow(graph))
  } else if (metric == 'meanMeanEdges') {
    return(mean(rowMeans(graph)))
  } else if (metric == 'meanDegree') {
    return(mean(rowSums(graph)))
  } else if (metric == 'medianMeanEdges') {
    return(median(rowMeans(graph)))
  } else if (metric == 'medianDegree') {
    return(median(rowSums(graph)))
  } else if (metric == 'meanShortestPath') {
    g <- graph_from_adjacency_matrix(graph, mode = 'undirected', weighted = NULL, diag = F)
    return(mean_distance(g, unconnected=T))
  } else if (metric == 'propDeg1') {
    return(mean(rowSums(graph) == 1))
  }
  # else if (metric == 'modularity') {
  #   return(modularity(graph_from_adjacency_matrix(graph, mode ='undirected', diag=F)))
  # }
}

#' Construct network using p-value and correlation thresholding and PNS (Alg 1)
#' 
#' @param X Gene expression matrix, samples x genes
#' @param pvals Vector of p-values
#' @param pvalthresh Gene p-value threshold
#' @param corthresh Neighbor correlation threshold
#' @param pvalscreen Vector of gene names corresponding to p-value thresholding (screen1_genes), default=NULL
#' @param lambdas Vector of regularization parameters
#' @param metrics Vector of metrics to calculate, default=NULL
#' @param degcut Drop nodes with degree less than degcut, default=2
#' @param maxnumcomp Maximum number of connected components (drop smallest components until you reach this amount), default=NULL
#' @param mincompsize Minimum component size (drop components smaller than this), default=5
#' @param degdivide Divide igraph node sizes by degdivide, default=10
#' @param standardize Whether or not to standardize data for PNS, default=T
#' @param savefile String denoting path and file name for where to store file, '_lambda<lambda>.rds' appended at end, default=NULL
#' @param saveviz String denoting path and file name for where to store graph viz, '<i>_cleangraph_lambda<lambda>.png' appended at end, default=NULL
#' @param verbose Whether or not to print out progress messages, default=T
#' @returns List of metric values (for each metric in metrics) for both pns and clean graphs for each lambda in lambdas
construct_graphs_metrics <- function(X, pvals, pvalthresh, corthresh, pvalscreen, 
                                     lambdas, metrics=NULL, 
                                     degcut=1, maxnumcomp=NULL, mincompsize=10, degdivide=10,
                                     standardize=T, savefile=NULL, saveviz=NULL, verbose=T) {
  
  # P-val screening and neighbor retrieval (Alg 1, Steps 1-3)
  pns_inputs <- screen_retrieve(X, pvals, pvalthresh, corthresh, pvalscreen, verbose=T)
  
  if (is.null(metrics)) {
    metrics <- c('R2', 'numEdges', 'numNodes', 'meanMeanEdges', 'meanDegree', 'medianMeanEdges', 'medianDegree', 'meanShortestPath', 'propDeg1')
  }
  if (any(!(metrics %in% c('R2', 'numEdges', 'numNodes', 'meanMeanEdges', 'meanDegree', 'medianMeanEdges', 'medianDegree', 'meanShortestPath', 'propDeg1')))) {
    stop("At least one metric is invalid, must only be c('R2', 'numEdges', 'numNodes', 'meanMeanEdges', 'meanDegree', 'medianMeanEdges', 'medianDegree', 'meanShortestPath', 'propDeg1').")
  }
  
  metric_list <- list()
  for (i in 1:length(lambdas)) {
    lambda <- lambdas[i]
    cat('lambda=', lambda, '\n')
    if (is.null(savefile)) {
      savefile_pnslambda <- NULL
      savefile_cleanlambda <- NULL
    } else {
      savefile_pnslambda <- paste0(savefile, 'pns_lambda', lambda, '.rds')
      savefile_cleanlambda <- paste0(savefile, 'clean_lambda', lambda, '.rds')
    }
    
    # Construct graph via PNS (Alg 1, Step 4)
    pnsgraph <- tryCatch(
      {
        run_pns(pns_inputs$X, pns_inputs$core_genes, pns_inputs$nbr_genes, lambda, 
                standardize=standardize, savefile=savefile_pnslambda, verbose=verbose)        
      }, error = function(msg){
        return(NULL)
      })
    
    cleangraph <- tryCatch(
      {
        clean_graph(pnsgraph$graphfinal, degcut=degcut, maxnumcomp=maxnumcomp, mincompsize=mincompsize, degdivide=degdivide, savefile=savefile_cleanlambda, verbose=T)
      }, error = function(msg){
        return(NULL)
      })
    
    if (all(rowSums(pnsgraph$graphfinal) == 0)) {
      pnsgraph <- NULL
      cleangraph <- NULL
    } 
    
    if (is.null(pnsgraph)) {
      pnsnet <- NULL
    } else {
      pnsnet <- pnsgraph$graphfinal
      # Remove lone floating nodes
      pnsnet <- pnsnet[which(rowSums(pnsnet) != 0), which(rowSums(pnsnet) != 0)]
    }
    
    if (is.null(cleangraph)) {
      cleannet <- NULL
    } else {
      cleannet <- cleangraph$ne
      
      if (!is.null(saveviz)) {
        saveviz_lambda <- paste0(saveviz, i,  '_clean_lambda', lambda, '.png')
        set.seed(42)
        if (length(V(cleangraph$net)) > 1750) {
          layout <- layout_with_drl(cleangraph$net)
        } else {
          layout <- layout_with_fr(cleangraph$net)
        }
        png(saveviz_lambda, width=2400, height=2400, res=300)
        par(mar=c(0, 0, 0, 0)+.1)
        plot.igraph(cleangraph$net,
                    layout = layout,
                    vertex.label=NA)
        if (!is.null(dev.list())) {
          dev.off()
        }
      }
    }
    
    for (metric in metrics) {
      if (is.null(pnsnet)) {
        pns_val <- NA
      } else {
        pns_val <-  tryCatch(
          {
            calc_metric(pnsnet, metric)          
          }, error = function(msg){
            return(NA)
          })
      }
      metric_list[[paste0('pns_', metric)]] <- c(metric_list[[paste0('pns_', metric)]], pns_val)
      
      if (is.null(cleannet)) {
        clean_val <- NA
      } else {
        clean_val <-  tryCatch(
          {
            calc_metric(cleannet, metric)          
            }, error = function(msg){
            return(NA)
          })
      }
      metric_list[[paste0('clean_', metric)]] <- c(metric_list[[paste0('clean_', metric)]], clean_val)
    }
  }
  return(metric_list)
}

#' #' Function to run PNS (Alg 1, Step 4) with K-fold CV for mean 
#' #' 
#' #' @param X Gene expression matrix of genes in set V, samples x genes
#' #' @param core_genes Vector containing the genes in set S
#' #' @param nbr_genes Vector containing the genes in set V / S
#' #' @param lambdas Vector of regularization parameters
#' #' @param standardize Whether or not to standardize data, default=T
#' #' @param savefile String denoting full path and file name for where to store file (rds), "_lambdaL" will be added, default=NULL
#' #' @param verbose Whether or not to print out progress messages, default=T
#' #' @returns A list containing
#' #' - genesfinal: Vector of final genes in network
#' #' - graphfinal: 0-1 Matrix denoting final gene network
#' #' - coef_graphfinal: Weighted matrix with coef values as edges 
#' run_pns_lambdas <- function(X, core_genes, nbr_genes, lambdas, K=5, standardize=T, savefile=NULL, verbose=F) {
#'   
#'   numcg <- length(core_genes)
#'   numng <- length(nbr_genes)
#'   
#'   if (numcg == 0) {
#'     stop('No core genes inputted (`core_genes` cannot be empty).')
#'   }
#'   if (numng == 0) {
#'     warning('No neighbor genes inputted (`nbr_genes` is empty).')
#'   }
#'   if (length(intersect(core_genes, nbr_genes))) {
#'     stop('`core_genes` and `nbr_genes` should have no intersection.')
#'   }
#'   if (!all(c(core_genes, nbr_genes) %in% colnames(X))) {
#'     stop('Not all `core_genes` and `nbr_genes` are in `X`.')
#'   }
#'   
#'   genes <- c(core_genes, nbr_genes)
#'   # Expression data of core genes
#'   yexp <- X[, core_genes] 
#'   # Expression data of neighbor genes 
#'   if (numng > 0) {
#'     xexp <- X[, nbr_genes] 
#'     xyexp <- cbind(yexp, xexp)
#'   } else{
#'     xyexp <- yexp
#'   }
#'   
#'   # Store mean cross-validated errors
#'   cvms_mx <- matrix(NA, nrow=length(core_genes), ncol=length(lambdas))
#'   for (i in 1:length(core_genes)) {
#'     cv_model <- cv.glmnet(x=xyexp[,-i], 
#'                           y=yexp[,i], 
#'                           family = "gaussian", lambda=lambdas, 
#'                           standardize = standardize, nfolds=K)
#'     cvms_mx[i,] <- cv_model$cvm
#'   }
#'   cvms <- colMeans(cvms_mx)
#'   
#'   if (!is.null(savefile)) {
#'     savepath <- strsplit(savefile, '.rds')
#'     
#'     for (lambda in lambdas) {
#'       if (verbose) {
#'         cat(paste0('lambda=', lambda, '\n'))
#'       }
#'       edge_coef <- function(i) {
#'         vectorr <- rep(0, numcg + numng)
#'         model <- glmnet(xyexp[,-i], 
#'                         yexp[,i], 
#'                         family = "gaussian", lambda = lambda, 
#'                         standardize = standardize)
#'         vectorr[-i] <- c(model$beta[, 1])
#'         return(vectorr)
#'       }
#'       
#'       edgecoef_mx <- sapply(1:numcg, edge_coef)
#'       if (numng > 0) {
#'         edgecoef_mx <- cbind(edgecoef_mx, matrix(0, nrow=numcg+numng, ncol=numng))
#'       }
#'       rownames(edgecoef_mx) = colnames(edgecoef_mx) = genes
#'       
#'       # Binarize matrix - E matrix
#'       binary_mx <- edgecoef_mx
#'       binary_mx <- (abs(binary_mx) > 0) * 1
#'       # Symmetrize matrix - Omega matrix
#'       binary_mx <- 1 - (1 - binary_mx) * (1 - t(binary_mx))
#'       
#'       final_genes <- genes #genes[which(rowSums(binary_mx) != 0)]
#'       final_mx <- binary_mx[final_genes, final_genes]
#'       finalcoef_mx <- edgecoef_mx[final_genes, final_genes]
#'       
#'       if (verbose) {
#'         cat(paste0('Number of edges among core genes: ', sum(final_mx[core_genes, core_genes])/2, '\n'))
#'       }
#'       
#'       final <- list()
#'       final$genesfinal <- final_genes
#'       final$graphfinal <- final_mx
#'       final$coef_graphfinal <- finalcoef_mx
#'       
#'       # Save final object
#'       saveRDS(final, paste0(savepath, '_lambda', lambda, '.rds'))
#'     }
#'   }
#'   
#'   return(cvms)
#' }

#' Function to clean/prune graph
#' 
#' @param graph 0-1 Adjacency matrix 
#' @param degcut Drop nodes with degree less than degcut, default=2
#' @param maxnumcomp Maximum number of connected components (drop smallest components until you reach this amount), default=NULL
#' @param mincompsize Minimum component size (drop components smaller than this), default=5
#' @param degdivide Divide igraph node sizes by degdivide, default=10
#' @param savefile String denoting full path and file name for where to store file (rds), default=NULL
#' @param verbose Whether or not to print out progress messages, default=T
#' @returns A list containing
#' - ne: Pruned 0-1 ajacency matrix
#' - net: Pruned igraph object
clean_graph <- function(graph, degcut=2, maxnumcomp=NULL, mincompsize=5, degdivide=10, savefile=NULL, verbose=T){
  
  if (nrow(graph) != ncol(graph)) {
    stop('`graph` must be square matrix.')
  }
  if (!isSymmetric(graph)) {
    stop('`graph` must be symmetric.')
  }
  if ((degcut <= 0) | (degdivide <= 0)) {
    stop('`degcut` must be > 0 and `degdivide` must be > 0.')
  }
  if (!is.null(savefile)) {
    savefilesplit <- strsplit(savefile, '\\.')[[1]]
    if (tail(savefilesplit, 1) != 'rds') {
      warning('`savefile` does not end in .rds, manually appending .rds.')
      savefile <- paste0(savefile, '.rds')
    }
  }
    
  if (verbose) {
    cat(paste0('Starting number of nodes: ', dim(graph)[1], '\n'))
  }
  net <- graph_from_adjacency_matrix(graph, mode = "undirected")
  net <- igraph::simplify(net, remove.multiple = T, remove.loops = T) 
  if (verbose) {
    cat(paste0('Number of starting solo island nodes: ', sum(degree(net) == 0), '\n'))
  }
  num_nodes_1 <- length(V(net))
  # Drop nodes with degree < `degcut`
  net <- delete_vertices(net, igraph::degree(net) < degcut)
  num_nodes_2 <- length(V(net))
  
  if (verbose) {
    cat(paste0('Dropping nodes with degree less than `degcut`=', degcut, '...'), '\n')
    cat(paste0('Number of nodes dropped: ', num_nodes_1 - num_nodes_2, '\n'))
    cat(paste0('Number of resulting solo island nodes (dropped): ', sum(degree(net) == 0), '\n'))
  }
  
  # Drop solo island nodes (nodes with no edges)
  net <- delete_vertices(net, igraph::degree(net) == 0)
  
  # Drop components with size smaller than `mincompsize`
  if (!is.null(mincompsize)) {
    comps <- components(net)
    dropcomp_idxs <- which(comps$csize < mincompsize)
    dropcomp_nodes <- names(comps$membership[comps$membership %in% dropcomp_idxs])
    net <- delete_vertices(net, dropcomp_nodes)
    if (verbose) {
      cat(paste0('Dropping components with size smaller than `mincompsize`=', mincompsize, '...'), '\n')
      cat(paste0('Number of components dropped: ', length(dropcomp_idxs)), '\n')
      cat(paste0('Number of nodes dropped: ', length(dropcomp_nodes)), '\n')
      cat(paste0('Number of resulting components: ', length(comps$csize) - length(dropcomp_idxs)), '\n')
    }
  }
  
  # Drop components to reach `maxnumcomp` components
  if (!is.null(maxnumcomp)) {
    comps <- components(net)
    numcomps <- length(comps$csize)
    if (numcomps > maxnumcomp) {
      numtodrop <- numcomps-maxnumcomp
      dropcomp_idxs <- order(comps$csize)[1:numtodrop]
      dropcomp_nodes <- names(comps$membership[comps$membership %in% dropcomp_idxs])
      net <- delete_vertices(net, dropcomp_nodes)
    }
    if (verbose) {
      cat(paste0('Dropping components to reach `maxnumcomp`=', maxnumcomp,' components...'), '\n')
      cat(paste0('Number of components dropped: ', numtodrop), '\n')
      cat(paste0('Number of nodes dropped: ', length(dropcomp_nodes)), '\n')
    }
  }
  
  if (verbose) {
    cat(paste0('Final number of nodes: ', length(V(net))), '\n')
    cat(paste0('Final number of edges: ', length(E(net))), '\n')
  }
  
  # Set node sizes (for plotting)
  deg <- igraph::degree(net, mode="all")
  V(net)$size <- deg/degdivide
  nodes <- V(net)$name
  
  res <- list()
  res$ne <- graph[nodes, nodes]
  res$net <- net
  
  # Save final object
  if (!is.null(savefile)) {
    saveRDS(res, savefile)
  }
  return(res)
}

#' Compute coefficient of determination (R^2) = square of the correlation coefficient for graph
#' How well does the degree distribution follow a power-law?
#' 
#' @param graph 0-1 symmetric adjacency matrix (output from construct_graph and/or clean_graph)
#' @returns Squared correlation value for input graph
get_R2 <- function(graph) {
  
  if (!isSymmetric(graph)) {
    stop('`graph` must be symmtric.')
  }
  
  pk <- table(colSums(graph))
  if (names(pk)[1]=='0') {
    pk <- pk[-1]
  }
  stopifnot(length(pk) > 0)
  k <- as.numeric(names(pk))
  R2 <- cor(log(pk), log(k))^2
  return(R2)
}
