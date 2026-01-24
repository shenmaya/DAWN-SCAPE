# Functions that should be in this file:
# hmrf main function
# pseudo conditional likelihood - PLCK

#' Calculate the log pseudo conditional likelihood. 
#' log prod_i p(I_i | {I_j}_{j in N_i}) = sum_i log(p(I_i | {I_j}_{j in N_i}))
#' 
#' @param genes Vector of gene names
#' @param dezscores Vector of DE z-scores, if not named, assumed to be in same order as `genes`
#' @param genzscores Vector of genetic z-scores, if not named, assumed to be in same order as `genes`
#' @param graph 0-1 matrix representing the graph
#' @param deinit DE z-score cutoff for initializing states (I)
#' @param geninit Genetic z-score cutoff for initializing states (I)
#' @param numiter Number of iterations to run 
#' @param b0_lb Lower bound for b0 parameters, must be negative (<0), default=-20
#' @param b1_ub Upper bound for b1 parameters, must be positive (>0), default=10
#' @param savefile String denoting full path and file name for where to store file (rds), default=NULL
#' @param verbose Whether or not to print out progress messages, default=T
#' @returns ...
#' @examples
#' 
run_hmrf <- function(genes, dezscores, genzscores, graph, deinit, geninit, detrim, gentrim,
                                                numiter, b0_lb=-20, b1_ub=10, savefile=NULL, verbose=T) {
  
  # Checks
  genes_len <- length(genes)
  dez_len <- length(dezscores)
  genz_len <- length(genzscores)
  graph_nrow <- nrow(graph)
  graph_ncol <- ncol(graph)
  if (length(unique(c(genes_len, dez_len, genz_len, graph_nrow, graph_ncol))) != 1) {
    stop('Check dimensions of inputs. `genes`, `dezscores`, `genzscores` (length), and `graph` (nrow, ncol) should have the sane dims.')
  }
  if (any(is.na(c(dezscores, genzscores, graph, deinit, geninit)))) {
    stop('Check inputs for NAs.')
  }
  if (any(is.infinite(c(dezscores, genzscores, deinit, geninit)))) {
    stop('Check inputs for Inf values.')
  }
  if (any(is.infinite(graph))) {
    stop('Check `graph` for Inf values.')
  }
  if ((b0_lb >= 0) | (b1_ub <= 0)) {
    stop('b0_lb must be less than 0 and b1_ub must be greater than 0.')
  }
  if (numiter <= 0) {
    stop('Invalid input of `numiter`, must be > 0, decimals will be floor-ed to nearest integer.')
  }
  if (!is.null(names(dezscores)) & !is.null(names(genzscores))) {
    if (length(intersect(genes, intersect(names(dezscores), names(genzscores)))) != genes_len) {
      stop('Check names of `dezscores` and `genzscores`, they do not match with `genes`.')
    } else {
      # Order z-scores by `genes`
      dezscores <- dezscores[genes]
      genzscores <- genzscores[genes]
    }
  }
  if (!is.null(savefile)) {
    savefilesplit <- strsplit(savefile, '\\.')[[1]]
    if (tail(savefilesplit, 1) != 'rds') {
      warning('`savefile` does not end in .rds, manually appending .rds.')
      savefile <- paste0(savefile, '.rds')
    }
    if (file.exists(savefile)) {
      cat('Warning: `savefile` already exists and will be overwritten. \n')
    }
  }
  
  
  if (verbose) {
    cat('Start HMRF analysis... \n')
  }
  
  # Initialize hidden state vector (Alg 2, Line 7, Eqn 7)
  Iupdate <- 1*((dezscores > deinit) & (genzscores > geninit)) + 2*((dezscores > deinit) & (genzscores <= geninit)) + 3*((dezscores <= deinit) & (genzscores > geninit)) + 4*((dezscores <= deinit) & (genzscores <= geninit))
  if (verbose) {
    print(table(Iupdate))
  }
  # old code - maybe we still do need this... especially for AD
  if (is.na(detrim)) {
    detrim <- Inf
  }
  if (is.na(gentrim)) {
    gentrim <- Inf
  }
  fixindex <- 1*((dezscores > detrim) & (genzscores > gentrim)) + 2*((dezscores > detrim) & (genzscores <= gentrim))+ 3*((dezscores <= detrim) & (genzscores > gentrim)) + 4*((dezscores <= detrim) & (genzscores <= gentrim))
  # fixindex <- 1*((dezscores > detrim) & (genzscores > gentrim)) + 2*((dezscores > detrim) & (genzscores <= gentrim))+ 3*(dezscores <= detrim) 
  print(table(fixindex))
  # fixindex <- rep(3, length(genes))
  
  # Initialize hidden variable-related parameters (Alg 2, Line 7, Eqn 8)
  b01 <- b02 <- b03 <- b11 <- b12 <- 0
  
  # Initialize normal parameters (Alg 2, Line 7, Eqn 8)
  mu1 <- mean(dezscores[(Iupdate %in% c(1, 2)) & (fixindex==4)])
  sigma1 <- mean((dezscores[Iupdate %in% c(1, 2)] - mu1)^2) # sd(dezscores[Iupdate %in% c(1, 2)])^2
  
  mu2 <- mean(genzscores[(Iupdate %in% c(1, 3)) & (fixindex==4)])
  sigma2 <- mean((genzscores[Iupdate %in% c(1, 3)] - mu2)^2) # sd(genzscores[Iupdate %in% c(1, 3)])^2
  
  sigma01 <- mean(dezscores[Iupdate %in% c(3, 4)]^2) # sd(dezscores[Iupdate %in% c(3, 4)])^2
  sigma02 <- mean(genzscores[Iupdate %in% c(2, 4)]^2) # sd(genzscores[Iupdate %in% c(2, 4)])^2
  
  # Clamping to try to prevent sigmas from becoming too small --> NAs and errored runs
  sigma1 <- max(sigma1, 1e-6)
  sigma2 <- max(sigma2, 1e-6)
  sigma01 <- max(sigma01, 1e-6)
  sigma02 <- max(sigma02, 1e-6)
  
  if (verbose) {
    cat("iteration", 0, ": b01 =", b01, ", b02 =", b02, ", b03 =", b03, ", b11 =", b11, ", b12 =", b12,
        ", mu1 =", mu1, ", sigma1 =", sigma1, ", mu2 =", mu2, ", sigma2 =", sigma2, ", sigma01 =", sigma01, ", sigma02 =", sigma02, "\n")
  }
  param_hist <- c(b01, b02, b03, b11, b12, mu1, sigma1, mu2, sigma2, sigma01, sigma02)

  state_hist <- Iupdate
  post1_hist <- c()
  post2_hist <- c()
  post3_hist <- c()
  post4_hist <- c()
  logppost_hist <- c()
  logpclk_hist <- c()
  b_osc_bools <- rep(F, 11)
  for (iter in 1:numiter) {
    
    # Save current hidden state vector and parameters
    Ibefore <- Iupdate
    b01_b <- b01 
    b02_b <- b02
    b03_b <- b03
    b11_b <- b11
    b12_b <- b12
    mu1_b <- mu1
    sigma1_b <- sigma1
    mu2_b <- mu2
    sigma2_b <- sigma2
    sigma01_b <- sigma01
    sigma02_b <- sigma02
    
    # Update hidden variable-related parameters (Alg 2, Line 10)
    numactive <- as.vector((Ibefore==1) %*% graph)
    numreactive <- as.vector((Ibefore==2) %*% graph)
    
    # Optimization code, switches to basic `maximize_logpclk` in case of error in `optim`
    b_bs <- c(b01_b, b02_b, b03_b, b11_b, b12_b)
    optim_result <- try({
      optim(
        par = b_bs,
        fn = logpclk_optim,
        states = Ibefore,
        numactive = numactive,
        numreactive = numreactive,
        method = "L-BFGS-B",  # Bounded optimization
        lower = c(b0_lb, b0_lb, b0_lb, 0, 0),  # Lower bounds
        upper = c(0, 0, 0, b1_ub, b1_ub)       # Upper bounds
      )
    }, silent = TRUE)
    
    optim_failed <- FALSE
    if (inherits(optim_result, "try-error")) {
      optim_failed <- TRUE
    } else {
      if (is.null(optim_result$par) || !all(is.finite(optim_result$par))) {
        optim_failed <- TRUE
      }
    }
    
    if (!optim_failed) {
      # optim succeeded
      b01 <- optim_result$par[1]
      b02 <- optim_result$par[2]
      b03 <- optim_result$par[3]
      b11 <- optim_result$par[4]
      b12 <- optim_result$par[5]
    } else {
      # Fall back to coordinate ascent
      warning("optim() failed or returned invalid values. Falling back to maximize_logpclk().")
      
      b_params_values <- maximize_logpclk(
        Ibefore, numactive, numreactive,
        b01_b, b02_b, b03_b, b11_b, b12_b,
        b0_lb = b0_lb, b1_ub = b1_ub,
        numiter = 20
      )
      
      b01 <- b_params_values[1]
      b02 <- b_params_values[2]
      b03 <- b_params_values[3]
      b11 <- b_params_values[4]
      b12 <- b_params_values[5]
    }

    posterior1s <- rep(NA, length(genes))
    posterior2s <- rep(NA, length(genes))
    posterior3s <- rep(NA, length(genes))
    posterior4s <- rep(NA, length(genes))
    
    # Update hidden state vector via ICM (Alg 2, Line 11)
    # Save the the posterior prob for each gene i to calculate the pseudo-posterior 
    # log(p(I | Z)) = log(\prod_i p(I_i | Z)) = \sum_i log(p(I_i | Z))
    post_is <- c()
    clk_is <- c()
    
    # for (node_i in sample(1:length(genes), replace=F)) {
    for (node_i in 1:length(genes)) {
      # Gene i's state
      state_i <- Ibefore[node_i]
      # Gene i's neighbors
      graph_i <- graph[node_i, ]
      # Number of active neighbors of gene i
      numactive_i <- as.numeric((Ibefore==1) %*% graph_i)
      # Number of reactive neighbors of gene i
      numreactive_i <- as.numeric((Ibefore==2) %*% graph_i)
      
      # Calculate CLK of each possible state (1, 2, 3, 4)
      qprobs <- clk(state_i, numactive_i, numreactive_i, b01, b02, b03, b11, b12, logbool=F, allbool=T)
      q1prob <- qprobs[1]
      q2prob <- qprobs[2]
      q3prob <- qprobs[3]
      q4prob <- qprobs[4]
      
      # Calculate the posterior probability of each possible state (1, 2, 3, 4)
      post1 <- dnorm(dezscores[node_i], mean=mu1, sd=sqrt(sigma1))*dnorm(genzscores[node_i], mean=mu2, sd=sqrt(sigma2))*q1prob
      post2 <- dnorm(dezscores[node_i], mean=mu1, sd=sqrt(sigma1))*dnorm(genzscores[node_i], mean=0, sd=sqrt(sigma02))*q2prob
      post3 <- dnorm(dezscores[node_i], mean=0, sd=sqrt(sigma01))*dnorm(genzscores[node_i], mean=mu2, sd=sqrt(sigma2))*q3prob
      post4 <- dnorm(dezscores[node_i], mean=0, sd=sqrt(sigma01))*dnorm(genzscores[node_i], mean=0, sd=sqrt(sigma02))*q4prob
      postsum <- post1+post2+post3+post4
      
      posterior1s[node_i] <- post1/postsum
      posterior2s[node_i] <- post2/postsum
      posterior3s[node_i] <- post3/postsum
      posterior4s[node_i] <- post4/postsum
      
      # Update with state with the largest posterior probability
      Iupdate[node_i] <- which.max(c(posterior1s[node_i], posterior2s[node_i], posterior3s[node_i], posterior4s[node_i]))
      
      # Save pseudo-posterior of gene i
      post_i <- c(posterior1s[node_i], posterior2s[node_i], posterior3s[node_i], posterior4s[node_i])[Iupdate[node_i]]
      post_is <- c(post_is, post_i)
      # Save pseudo-conditional likelihood of gene i
      clk_i <- c(q1prob, q2prob, q3prob, q4prob)[Iupdate[node_i]]
      clk_is <- c(clk_is, clk_i)
    }
    logppost <- sum(log(post_is))
    logpclk <- sum(log(clk_is))
    
    # Update normal parameters (Alg 2, Line 12)
    mu1 <- sum((dezscores*(posterior1s+posterior2s))[fixindex==4])/sum((posterior1s+posterior2s)[fixindex==4])
    sigma1 <- sum((((dezscores - mu1)^2)*(posterior1s+posterior2s))[fixindex==4])/sum((posterior1s+posterior2s)[fixindex==4])
    
    mu2 <- sum((genzscores*(posterior1s+posterior3s))[fixindex==4])/sum((posterior1s+posterior3s)[fixindex==4])
    sigma2 <- sum((((genzscores - mu2)^2)*(posterior1s+posterior3s))[fixindex==4])/sum((posterior1s+posterior3s)[fixindex==4])
    
    sigma01 <- sum(((dezscores^2)*(posterior3s+posterior4s))[fixindex==4])/sum((posterior3s+posterior4s)[fixindex==4])
    sigma02 <- sum(((genzscores^2)*(posterior2s+posterior4s))[fixindex==4])/sum((posterior2s+posterior4s)[fixindex==4])
    
    # Clamping to try to prevent sigmas from becoming too small --> NAs and errored runs
    sigma1 <- max(sigma1, 1e-6)
    sigma2 <- max(sigma2, 1e-6)
    sigma01 <- max(sigma01, 1e-6)
    sigma02 <- max(sigma02, 1e-6)
    
    # # Check for oscillation in b params
    # # values which should perhaps be included as arg: window, tol, threshold, alpha
    # param_hist_temp <- rbind(param_hist, c(b01, b02, b03, b11, b12, mu1, sigma1, mu2, sigma2, sigma01, sigma02))
    # b_alphas <- rep(0.8, 11)
    # # Continue dampening flagged oscillating b parameters and see if there are any new ones
    # b_osc_bools_curr <- apply(param_hist_temp, 2, detect_oscillation, window=10)
    # b_osc_bools <- (b_osc_bools | b_osc_bools_curr)
    # # If b parameter not oscillating, do not dampen
    # b_alphas[!b_osc_bools] <- 0
    # bs_curr <- c(b01, b02, b03, b11, b12, mu1, sigma1, mu2, sigma2, sigma01, sigma02) # c(b01, b02, b03, b11, b12)
    # bs_prev <- c(b01_b, b02_b, b03_b, b11_b, b12_b, mu1_b, sigma1_b, mu2_b, sigma2_b, sigma01_b, sigma02_b) # c(b01_b, b02_b, b03_b, b11_b, b12_b)
    # bs_dampened <- (bs_curr * (1-b_alphas)) + (bs_prev * b_alphas)
    # 
    # param_hist <- rbind(param_hist, c(bs_dampened[1], bs_dampened[2], bs_dampened[3], bs_dampened[4], bs_dampened[5], 
    #                                   bs_dampened[6], bs_dampened[7], bs_dampened[8], bs_dampened[9], 
    #                                   bs_dampened[10], bs_dampened[11]))
    
    param_hist <- rbind(param_hist, c(b01, b02, b03, b11, b12, mu1, sigma1, mu2, sigma2, sigma01, sigma02))
    if (verbose) {
      cat("iteration", iter, ": b01 =", b01, ", b02 =", b02, ", b03 =", b03, ", b11 =", b11, ", b12 =", b12,
          ", mu1 =", mu1, ", sigma1 =", sigma1, ", mu2 =", mu2, ", sigma2 =", sigma2, ", sigma01 =", sigma01, ", sigma02 =", sigma02, "\n")
    }
    state_hist <- rbind(state_hist, Iupdate)
    post1_hist <- rbind(post1_hist, posterior1s)
    post2_hist <- rbind(post2_hist, posterior2s)
    post3_hist <- rbind(post3_hist, posterior3s)
    post4_hist <- rbind(post4_hist, posterior4s)
    logppost_hist <- c(logppost_hist, logppost)
    logpclk_hist <- c(logpclk_hist, logpclk)
  }
  
  colnames(param_hist) <- c('b01', 'b02', 'b03', 'b11', 'b12', 'mu1', 'sigma1', 'mu2', 'sigma2', 'sigma01', 'sigma02')
  rownames(param_hist) <- NULL
  rownames(state_hist) <- NULL
  colnames(post1_hist) <- genes
  rownames(post1_hist) <- NULL
  colnames(post2_hist) <- genes
  rownames(post2_hist) <- NULL
  colnames(post3_hist) <- genes
  rownames(post3_hist) <- NULL
  colnames(post4_hist) <- genes
  rownames(post4_hist) <- NULL
  
  res <- list()
  res$Iupdate <- Iupdate
  res$state_hist <- state_hist
  res$params <- list()
  res$params$b01 <- b01
  res$params$b02 <- b02
  res$params$b03 <- b03
  res$params$b11 <- b11
  res$params$b12 <- b12
  res$params$mu1 <- mu1
  res$params$sigma1 <- sigma1
  res$params$mu2 <- mu2
  res$params$sigma2 <- sigma2
  res$params$sigma01 <- sigma01
  res$params$sigma02 <- sigma02
  res$param_hist <- param_hist
  res$post1_hist <- post1_hist
  res$post2_hist <- post2_hist
  res$post3_hist <- post3_hist
  res$post4_hist <- post4_hist
  res$logppost_hist <- logppost_hist
  res$logpclk_hist <- logpclk_hist
  
  # Save final object
  if (!is.null(savefile)) {
    saveRDS(res, savefile)
  }
  return(res)
}

#' Calculate the log pseudo conditional likelihood. 
#' log prod_i p(I_i | {I_j}_{j in N_i}) = sum_i log(p(I_i | {I_j}_{j in N_i}))
#' 
#' @param args List containing 8 values:
#' - args[[1]]: initthres_de
#' - args[[2]]: initthres_gen
#' ...
#' @returns ...
#' @examples
#' 
run_hmrf_parallel <- function(args) {
  
  if (!exists(quote(zscorede))) {
    stop('zscorede variable not found. should be named array of DE z-scores.')
  } else if (!exists(quote(zscoregen))) {
    stop('zscoregen variable not found. should be named array of genetic z-scores.')
  } else if (!exists(quote(graph))) {
    stop('graph variable not found. should be igraph object. ')
  } else if (!exists(quote(layout))) {
    stop('layout variable not found. should be layout generated from igraph. ')
  }
  
  if (length(zscorede) != length(zscoregen)) {
    stop('zscorede and zscoregen should have the same length.')
  }
  if (nrow(layout) != length(zscorede)) {
    stop('layout should have a row per gene.')
  }
  
  initthres_de <- args[[1]]
  initthres_gen <- args[[2]]
  trimthres_de <- args[[3]] 
  trimthres_gen <- args[[4]] 
  numiter <- args[[5]]
  contrun <- args[[6]]
  savefile <- args[[7]]  
  saveviz <- args[[8]]
  
  if (tail(strsplit(savefile, '\\.')[[1]], 1) != 'rds') {
    warning('`savefile` does not end in .rds, manually appending .rds.')
    savefile <- paste0(savefile, '.rds')
  }
  if (length(args) > 8) {
    seed <- args[[9]]
    set.seed(seed)
    genelist <- sample(V(graph)$name, length(V(graph)$name), replace=FALSE)
    cat(paste0('seed ', seed, ':'), genelist[1:10], '\n')
  } else {
    genelist <- V(graph)$name
  }
  
  t0 <- Sys.time()
  hmrf_res <- run_hmrf(genes=genelist, dezscores=zscorede, genzscores=zscoregen, 
                       graph=as_adjacency_matrix(graph), 
                       deinit=initthres_de, geninit=initthres_gen, 
                       detrim=trimthres_de, gentrim=trimthres_gen,
                       numiter=numiter, savefile=savefile, verbose=T)
  t1 <- Sys.time()
  print(t1-t0)
  
  la <- hmrf_res$Iupdate
  la <- la[V(graph)$name]
  
  saveres <- readRDS(paste0(savefile))
  cols = c('#E66100','#1A85FF','#C5B6E0','#C5B6E0')
  pdf(saveviz, height = 7, width = 7)
  par(mar=c(0,0,0,0)+.1)
  plot.igraph(graph,
              layout=layout,
              vertex.size=V(graph)$size,
              vertex.color= cols[la],
              vertex.frame.color = cols[la],
              edge.width = 0.5,
              edge.color = 'gray',
              vertex.label=NA)
  title(paste0("state graph (iter: ", nrow(saveres$param_hist)-1, ") \n init_de: ", round(initthres_de, 3), ', init_gen: ', round(initthres_gen, 3),
               ', trimde: ', round(trimthres_de, 3), ', trimgen: ', round(trimthres_gen, 3)),line=-3)
  if (!is.null(dev.list())) {
    dev.off()
  }
}

#' Maximize the log pseudo conditional likelihood. 
#' argmax_{b0s, b1s} log prod_i p(I_i | {I_j}_{j in N_i}) = sum_i log(p(I_i | {I_j}_{j in N_i}))
#' 
#' @param states Vector of states (I)
#' @param numactive Vector of the number of active neighbors each gene has, assumed same order as `states`
#' @param numreactive Vector of the number of reactive neighbors each gene has, assumed same order as `states`
#' @param b01-b12 Current parameter estimates of b01, b02, b03, b11, b12
#' @param b0_lb Lower bound for b0 parameters, must be negative (<0), default=-20
#' @param b1_ub Upper bound for b1 parameters, must be positive (>0), default=10
#' @param numiter Number of iterations to run 
#' @returns ...
#' @examples
#' 
maximize_logpclk <- function(states, numactive, numreactive, 
                             b01, b02, b03, b11, b12, b0_lb=-20, b1_ub=10, numiter=20) {
  b01_iter <- b01
  b02_iter <- b02
  b03_iter <- b03
  b11_iter <- b11
  b12_iter <- b12
  
  for (k in 1:numiter) {
    b01_iter_k <- optimize(logpclk, c(b0_lb,0), b02 = b02_iter, b03 = b03_iter, 
                           b11 = b11_iter, b12 = b12_iter,
                           states = states, numactive = numactive, numreactive = numreactive,
                           maximum = T)$maximum
    b02_iter_k <- optimize(logpclk, c(b0_lb,0), b01 = b01_iter_k, b03 = b03_iter, 
                           b11 = b11_iter, b12 = b12_iter,
                           states = states, numactive = numactive, numreactive = numreactive,
                           maximum = T)$maximum
    b03_iter_k <- optimize(logpclk, c(b0_lb,0), b01 = b01_iter_k, b02 = b02_iter_k, 
                           b11 = b11_iter, b12 = b12_iter,
                           states = states, numactive = numactive, numreactive = numreactive,
                           maximum = T)$maximum
    b11_iter_k <- optimize(logpclk, c(0, b1_ub), b01 = b01_iter_k, b02 = b02_iter_k, b03 = b03_iter_k, 
                           b12 = b12_iter,
                           states = states, numactive = numactive, numreactive = numreactive,
                           maximum = T)$maximum
    b12_iter_k <- optimize(logpclk, c(0, b1_ub), b01 = b01_iter_k, b02 = b02_iter_k, b03 = b03_iter_k, 
                           b11 = b11_iter_k,
                           states = states, numactive = numactive, numreactive = numreactive,
                           maximum = T)$maximum
    if (abs(b01_iter_k - b01_iter) < 10^-5 & abs(b02_iter_k - b02_iter) < 10^-5 & abs(b03_iter_k - b03_iter) < 10^-5 & 
        abs(b11_iter_k - b11_iter) < 10^-5 & abs(b12_iter_k - b12_iter) < 10^-5) {
      break
    }
    b01_iter <- b01_iter_k
    b02_iter <- b02_iter_k
    b03_iter <- b03_iter_k
    b11_iter <- b11_iter_k
    b12_iter <- b12_iter_k
  }
  return(c(b01_iter, b02_iter, b03_iter, b11_iter, b12_iter))
}

# Define the function for optimization
logpclk_optim <- function(b_vec, states, numactive, numreactive) {
  # Extract parameters
  b01 <- b_vec[1]
  b02 <- b_vec[2]
  b03 <- b_vec[3]
  b11 <- b_vec[4]
  b12 <- b_vec[5]
  
  # Compute the negative log pseudo conditional likelihood (since optim minimizes)
  return(-logpclk(states, numactive, numreactive, b01, b02, b03, b11, b12))
}

# NOTES: is plck only used in log form? yes? 
#' Calculate the log pseudo conditional likelihood (CLK). 
#' log prod_i p(I_i | {I_j}_{j in N_i}) = sum_i log(p(I_i | {I_j}_{j in N_i}))
#' 
#' @param states Vector of states (I).
#' @param numactive Vector of the number of active neighbors each gene has, assumed same order as `states`.
#' @param numreactive Vector of the number of reactive neighbors each gene has, assumed same order as `states`.
#' @param b01-b12 Current parameter estimates of b01, b02, b03, b11, b12.
#' @param sumbool Whether the log CLKs should be summed, default=T.
#' @returns A scalar (sumbool=T) or vector (sumbool=F).
#' @examples
#' logpclk(c(1, 1, 3, 4, 2), c(11, 8, 2, 0, 2), c(3, 4, 5, 1, 9), -0.1, -0.2, -0.1, 2.1, 3.1)
logpclk <- function(states, numactive, numreactive, b01, b02, b03, b11, b12, sumbool=T) {
  clk_vals <- unname(mapply(clk, states, numactive, numreactive, 
                            MoreArgs = list(b01=b01, b02=b02, b03=b03, b11=b11, b12=b12,
                                            logbool=T, allbool=F)))
  if (sumbool) {
    return(sum(clk_vals))
  } else {
    return(clk_vals)
  }
}

#' Calculate the (log) conditional likelihood (CLK) of I_i given it's neighbors
#' 
#' @param state_i The state of node i, must be one of c(1, 2, 3, 4).
#' @param numactive_i The number of active neighbors of node i, integer.
#' @param numreactive_i The number of reactive neighbors of node i, integer.
#' @param b01-b12 Current parameter estimates of b01, b02, b03, b11, b12.
#' @param logbool Whether to return log CLK or CLK (boolean).
#' @param allbool Whether to return CLK for all states, default=F.
#' @returns A scalar (allbool=F) or vector (allbool=T).
#' @examples
#' clk(1, 10, 3, -0.1, -0.2, -0.1, 2.1, 3.1)
clk <- function(state_i, numactive_i, numreactive_i, b01, b02, b03, b11, b12, logbool, allbool=F) {
  
  q1 <- exp(b01 + (b11 * numactive_i))
  q2 <- exp(b02 + (b12 * numreactive_i))
  q3 <- exp(b03)
  q4 <- 1 # exp(0) = 1
  qsum <- q1 + q2 + q3 + q4
  qprobs <- c(q1, q2, q3, q4)/qsum
  
  if (logbool) {
    qprobs <- log(qprobs)
  }
  if (!allbool) {
    return(qprobs[state_i])
  } else {
    return(qprobs)
  }
}

#' Obtain final states by finding most likely state based on final parameters
#' 
#' @param genes Vector of gene names
#' @param graph 0-1 matrix representing the graph
#' @param Ibefore The last round of states (Iupdate from result of `run_hmrf`)
#' @param dezscores Vector of DE z-scores, if not named, assumed to be in same order as `genes`
#' @param genzscores Vector of genetic z-scores, if not named, assumed to be in same order as `genes`
#' @param b01-b12 Final parameters
#' @param mu1-sigma02 Final normal parameters
#' @returns A list containing: 
#' Iupdate - Vector of states (integers in c(1, 2, 3, 4))
#' logppost - Float of log pseudo-posterior (averaged over genes)
#' logpclk - Float of log conditional likelihood (averaged over genes)
#' @examples
final_states <- function(genes, graph, Ibefore, 
                         dezscores, genzscores,
                         b01, b02, b03, b11, b12, 
                         mu1, sigma1, mu2, sigma2,
                         sigma01, sigma02) {
  
  posterior1s <- rep(NA, length(genes))
  posterior2s <- rep(NA, length(genes))
  posterior3s <- rep(NA, length(genes))
  posterior4s <- rep(NA, length(genes))
  
  post_is <- c()
  clk_is <- c()
  
  Iupdate <- rep(0, length(genes))
  # for (node_i in sample(1:length(genes), replace=F)) {
  for (node_i in 1:length(genes)) {
    # Gene i's state
    state_i <- Ibefore[node_i]
    # Gene i's neighbors
    graph_i <- graph[node_i, ]
    # Number of active neighbors of gene i
    numactive_i <- as.numeric((Ibefore==1) %*% graph_i)
    # Number of reactive neighbors of gene i
    numreactive_i <- as.numeric((Ibefore==2) %*% graph_i)
    
    # Calculate CLK of each possible state (1, 2, 3, 4)
    qprobs <- clk(state_i, numactive_i, numreactive_i, b01, b02, b03, b11, b12, logbool=F, allbool=T)
    q1prob <- qprobs[1]
    q2prob <- qprobs[2]
    q3prob <- qprobs[3]
    q4prob <- qprobs[4]
    
    # Calculate the posterior probability of each possible state (1, 2, 3, 4)
    post1 <- dnorm(dezscores[node_i], mean=mu1, sd=sqrt(sigma1))*dnorm(genzscores[node_i], mean=mu2, sd=sqrt(sigma2))*q1prob
    post2 <- dnorm(dezscores[node_i], mean=mu1, sd=sqrt(sigma1))*dnorm(genzscores[node_i], mean=0, sd=sqrt(sigma02))*q2prob
    post3 <- dnorm(dezscores[node_i], mean=0, sd=sqrt(sigma01))*dnorm(genzscores[node_i], mean=mu2, sd=sqrt(sigma2))*q3prob
    post4 <- dnorm(dezscores[node_i], mean=0, sd=sqrt(sigma01))*dnorm(genzscores[node_i], mean=0, sd=sqrt(sigma02))*q4prob
    postsum <- post1+post2+post3+post4
    
    posterior1s[node_i] <- post1/postsum
    posterior2s[node_i] <- post2/postsum
    posterior3s[node_i] <- post3/postsum
    posterior4s[node_i] <- post4/postsum
    
    # Update with state with the largest posterior probability
    Iupdate[node_i] <- which.max(c(posterior1s[node_i], posterior2s[node_i], posterior3s[node_i], posterior4s[node_i]))
    
    # Save pseudo-posterior of gene i
    post_i <- c(posterior1s[node_i], posterior2s[node_i], posterior3s[node_i], posterior4s[node_i])[Iupdate[node_i]]
    post_is <- c(post_is, post_i)
    # Save pseudo-conditional likelihood of gene i
    clk_i <- c(q1prob, q2prob, q3prob, q4prob)[Iupdate[node_i]]
    clk_is <- c(clk_is, clk_i)
  }
  logppost <- sum(log(post_is))
  logpclk <- sum(log(clk_is))
  
  return(list(Iupdate=Iupdate,
              logppost=logppost,
              logpclk=logpclk))
}