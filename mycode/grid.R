library(ggplot2)
library(grid)
library(gridExtra)

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
  hmrf_res <- run_hmrf(genelist, zscorede, zscoregen, as_adjacency_matrix(graph), initthres_de, initthres_gen, 
                       numiter, b0_lb=-50, b1_ub=50, savefile=savefile, verbose=T)
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
run_hmrf_parallel_oldversion <- function(args) {
  
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
  hmrf_res <- run_hmrf_oldversion(genelist, zscorede, zscoregen, as_adjacency_matrix(graph), initthres_de, initthres_gen, 
                       numiter, savefile=savefile, verbose=T)
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
  title(paste0("state graph (iter: ", nrow(saveres$param_hist), ") \n init_de: ", round(initthres_de, 3), ', init_gen: ', round(initthres_gen, 3),
               ', trimde: ', round(trimthres_de, 3), ', trimgen: ', round(trimthres_gen, 3)),line=-3)
  if (!is.null(dev.list())) {
    dev.off()
  }
}

viz_mx_pval <- function(val_mx, pval_seq_de, pval_seq_gen, title, saveviz, maxdim=2000, fontsize=10, null_mx=NULL) {
  
  if (!dir.exists(dirname(saveviz))) {
    stop("invalid file path!")
  }
  if (!is.null(null_mx)) {
    if (!all(dim(val_mx) == dim(null_mx))) {
      stop("val_mx and null_mx need to have the same dimensions.")
    }
  }
  
  if (!is.null(null_mx)) {
    val_mx <- val_mx * null_mx
    
  }
  width <- maxdim*(length(pval_seq_de)/max(length(pval_seq_de), length(pval_seq_gen)))
  height <- maxdim*(length(pval_seq_gen)/max(length(pval_seq_de), length(pval_seq_gen)))
  
  jpeg(paste0(saveviz), width = width+50, height = height, res = 300, quality = 100)
  my_title <- textGrob(title, gp = gpar(fontsize = fontsize, fontface = "bold"))
  ht <- pheatmap::pheatmap(val_mx, 
                           cluster_rows=FALSE, cluster_cols=FALSE, 
                           display_numbers = TRUE, number_color = "black",
                           labels_col = colnames(val_mx), labels_row = rownames(val_mx)) 
  
  grid::grid.newpage()
  grid::grid.draw(ht$gtable)
  grid.arrange(grobs = list(my_title, ht[[4]]), heights = c(0.1, 1))
  dev.off()
}

# Function to create heatmap of number/proportion of starting active/state 1 genes 
viz_startmx_pval <- function(zscorede, zscoregen, pval_seq_de, pval_seq_gen, title, savepath, height=2000, state=1) {
  
  if (length(zscorede) != length(zscoregen)) {
    stop('zscorede and zscoregen should have the same lengths.')
  }
  
  if (length(grep('.', savepath, fixed=TRUE)) >= 1) {
    warning('savepath should not contain file extension, e.g. ".jpg"', immediate. = TRUE)
    Sys.sleep(5)
  }
  
  state_str <- c('Active', 'Reactive', 'State 3', 'State 4')[state]
  
  # how many state ?? genes do you start with?
  states_mx <- matrix(NA, nrow=length(pval_seq_gen), ncol=length(pval_seq_de))
  # i: genetic, j: DE
  for (idx_i in 1:length(pval_seq_gen)) {
    p_i <- pval_seq_gen[idx_i]
    z_i <- qnorm(1-p_i)
    
    for (idx_j in 1:length(pval_seq_de)) {
      p_j <- pval_seq_de[idx_j]
      z_j <- qnorm(1-p_j)
      
      if (state == 1) {
        num_state <- sum((zscorede > z_j) & (zscoregen > z_i))
      } else if (state == 2) {
        num_state <- sum((zscorede > z_j) & (zscoregen <= z_i))
      } else if (state == 3) {
        num_state <- sum((zscorede <= z_j) & (zscoregen > z_i))
      } else {
        num_state <- sum((zscorede <= z_j) & (zscoregen <= z_i))
      }
      
      states_mx[idx_i, idx_j] <- num_state
    }
  }
  rownames(states_mx) <- paste0('gen_p', pval_seq_gen)
  colnames(states_mx) <- paste0('de_p', pval_seq_de)

  jpeg(paste0(savepath, "_count.jpg"), width = height + (round(height/2000) * 200), height = height, res = 300, quality = 100)
  my_title <- textGrob(paste0(title, ': Number of ', state_str, ' Starting Genes'), gp = gpar(fontsize = 8, fontface = "bold"))
  ht <- pheatmap::pheatmap(round(states_mx), 
                           cluster_rows=FALSE, cluster_cols=FALSE, 
                           display_numbers = TRUE, number_color = "black")
  grid::grid.newpage()
  grid::grid.draw(ht$gtable)
  grid.arrange(grobs = list(my_title, ht[[4]]), heights = c(0.1, 1))
  dev.off()
  
  jpeg(paste0(savepath, "_prop.jpg"), width = height + (round(height/2000) * 200), height = height, res = 300, quality = 100)
  my_title <- textGrob(paste0(title, ': Proportion of ', state_str, ' Starting Genes'), gp = gpar(fontsize = 8, fontface = "bold"))
  ht <- pheatmap::pheatmap(round(states_mx/length(zscorede), 3), 
                           cluster_rows=FALSE, cluster_cols=FALSE, 
                           display_numbers = TRUE, number_color = "black")
  grid::grid.newpage()
  grid::grid.draw(ht$gtable)
  grid.arrange(grobs = list(my_title, ht[[4]]), heights = c(0.1, 1))
  dev.off()
}

viz_paramconv_iter <- function(params_list, title, savepath, numiters=NULL) {
  
  if (length(grep('.', savepath, fixed=TRUE)) >= 1) {
    warning('savepath should not contain file extension, e.g. ".jpg"', immediate. = TRUE)
    Sys.sleep(5)
  }
  if (is.null(numiters)) {
    numiters <- params_list[[1]][[5]]
  }
  
  train_df <- data.frame(b01=double(), b02=double(), b03=double(), b11=double(), b12=double(),
                         mu1=double(), mu2=double(), sigma1=double(), sigma2=double(), 
                         sigma01=double(), sigma01=double(), 
                         iter=double(), paramnum=double(), paramstr=character())
  for (i in 1:length(params_list)) {
    
    file_path <- paste0(params_list[[i]][[7]])
    if (file.exists(file_path)) {
      
      hmrf_res <- readRDS(file_path)
      paramnum <- as.numeric(strsplit(tail(strsplit(params_list[[i]][[7]], '_')[[1]], 1), '.rds')[[1]])
      
      train_i_df <- hmrf_res$param_hist
      train_i_df <- cbind(cbind(train_i_df, 1:nrow(train_i_df)), rep(paramnum, nrow(train_i_df)))
      colnames(train_i_df) <- c(colnames(train_i_df)[1:11], 'iter', 'paramnum')
      train_df <- rbind(train_df, train_i_df)
      
    }
  }
  train_df$paramnum <- factor(train_df$paramnum, unique(train_df$paramnum))

  params <- c('b01', 'b02', 'b03', 'b11', 'b12', 'mu1', 'mu2', 'sigma1', 'sigma2', 'sigma01', 'sigma02')
  for (i in 1:length(params)) {
    train_df_i <- train_df
    train_df_i <- train_df_i[train_df_i$iter <= numiters,]
    ggplot(train_df_i, aes(x=iter, y=train_df_i[,params[i]], color=paramnum)) +
      geom_line() +
      labs(title = paste0(params[i], ': ', title),
           x = "Iteration",
           y = params[i],
           color = "Parameter Set") +
      theme_minimal()
    
    ggsave(paste0(savepath, '_', params[i], '.jpg'),
           width=15,
           height=6)
  }
}