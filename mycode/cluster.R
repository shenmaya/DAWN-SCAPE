# Functions for clustering
library(igraph)
library(RColorBrewer)

#' Leiden clustering
#' 
#' @param 
#' @returns 
cluster_net <- function(graph, X, k=10, resol=1, folder='./', 
                        hcfilename=NA, savefigs=T, ext='') {
  # Get initial Leiden clustering, igraph function
  med <- 'Leiden'
  com <- cluster_leiden(graph,
                       objective_function = c("CPM", "modularity"),
                       resolution = resol,
                       beta = 0.01,
                       initial_membership = NULL,
                       n_iterations = 2,
                       vertex_weights = NULL)
  labels <- membership(com)
  
  mapping <- sort.int(table(labels),index.return = TRUE,decreasing = TRUE)$ix
  labels <- sapply(labels, function(x) match(x,mapping))
  eg = eigengraph(X[V(graph)$name,], labels)
  
  sizecat = cut(as.numeric(table(labels)),
                breaks=c(0,5,10,20,50,100,200,500,1000),
                labels=c('<5', '5~10', '10~20', '20~50','50~100','100~200','200~500','>500'))
  colo = colorRampPalette(c("white", "darkgreen"))(length(levels(sizecat)))
  colorrow = c('<5'=colo[1], '5~10'=colo[2], '10~20'=colo[3], '20~50'=colo[4],'50~100'=colo[5],'100~200'=colo[6],'200~500'=colo[7],'>500'=colo[8])
  annotation_row = data.frame(size = sizecat)
  rownames(annotation_row) = as.character(1:length(unique(labels)))
  
  rownames(eg[[1]]) = colnames(eg[[1]])=1:length(unique(labels))
  
  if (savefigs) {
    pheatmap::pheatmap(eg[[1]], cluster_rows = FALSE, cluster_cols = FALSE, 
                       filename = paste0(folder, 'eg', med, '_resol', resol, '_k', k, ext, '.pdf'),
                       width = 10, height = 10, fontsize = 15, 
                       breaks = seq(-1,1,length.out=100),
                       cellwidth = 400/length(unique(labels)),
                       cellheight = 400/length(unique(labels)),
                       annotation_colors = list(size=colorrow))
  }
  
  d <- as.dist((1-eg[[1]])/2)   
  h <- hclust(d)
  if (!is.na(hcfilename)) {
    saveRDS(h, paste0(hcfilename))
  } 
  # if (is.na(hcfilename)) {
  #   saveRDS(h, paste0(folder, 'hclust_', med, '_resol', resol, '_k', k, ext, '.rds'))
  # } else {
  #   saveRDS(h, paste0(hcfilename))
  # }
  # If number of clusters/communities > k, cut dendrogram to get to k clusters
  if (k<length(unique(labels))){
    clusters = cutree(h, k=k)
  } else{
    clusters = 1:length(unique(labels))
  }
  
  
  #print(table(clusters))
  biglabels <- sapply(labels, function(x)clusters[x])
  colorplate = c('forestgreen', brewer.pal(n = max(length(unique(biglabels))-2, 3), name = "Set3"), 'black')
  colorplate = colorplate[1:length(unique(biglabels))]
  names(labels) = V(graph)$name
  names(biglabels) = V(graph)$name
  
  
  eg = eigengraph(X[V(graph)$name,], biglabels)
  sizecat = cut(as.numeric(table(biglabels)),
                breaks=c(0,5,10,20,50,100,200,500,1000),
                labels=c('<5', '5~10', '10~20', '20~50','50~100','100~200','200~500','>500'))
  colo = colorRampPalette(c("white", "darkgreen"))(length(levels(sizecat)))
  colorrow = c('<5'=colo[1], '5~10'=colo[2], '10~20'=colo[3], '20~50'=colo[4],'50~100'=colo[5],'100~200'=colo[6],'200~500'=colo[7],'>500'=colo[8])
  annotation_row = data.frame(size = sizecat)
  rownames(annotation_row) = as.character(1:length(unique(biglabels)))
  rownames(eg[[1]]) = colnames(eg[[1]])=1:length(unique(biglabels))
  if (savefigs) {
    pheatmap::pheatmap(eg[[1]], cluster_rows = FALSE, cluster_cols = FALSE, 
                       filename = paste0(folder,'eg',med,'_resol', resol,'_k', k, ext, '_aftermerge.pdf'),
                       width = 10, height = 10, fontsize = 20, 
                       breaks = seq(-1,1,length.out=100),
                       cellwidth = 300/length(unique(biglabels)),
                       cellheight = 300/length(unique(biglabels)),
                       annotation_row = annotation_row,
                       annotation_colors = list(size=colorrow))
  }
  
  return(list(labels = biglabels, colors = colorplate, sublabels = labels))
  
}

run_cluster_net_reps <- function(param, cleanraw, bulk, folder, base_seed = 1234, reps = 10, hcfilename = NA, savefigs = FALSE) {
  # Extract resol and k from param
  resol <- param$resol
  k <- param$k
  param_id <- param$id  # Make sure this is added to param_grid beforehand
  
  # Store all results
  results_list <- vector("list", reps)
  # Create folder where labels will be saved
  subfolder <- paste0(folder, 'paramsid', param_id, '_resol', resol, '_k', k, '/')
  dir.create(subfolder)
  
  for (i in seq_len(reps)) {
    # Set a unique seed for reproducibility
    set.seed(base_seed + i + reps * param_id)
    
    # Run clustering
    result <- cluster_net(
      graph = cleanraw,
      X = bulk,
      k = k,
      resol = resol,
      folder = folder,
      hcfilename = NA,
      savefigs = savefigs
    )
    modua <- modularity(cleanraw, result$labels)
    modu <- modularity(cleanraw, result$sublabels)
    smclust <- sum(table(result$labels) < 10)
    
    saveRDS(result, paste0(subfolder, 'labels_resol', resol, '_k', k, '_rep', i, '.rds'))

    results_list[[i]] <- data.frame(
      resol = resol,
      k = k,
      run_id = i,
      modua = modua,
      modu = modu,
      smclust = smclust
    )
  }
  
  # Combine results into a data.frame
  do.call(rbind, results_list)
}

eigengraph<-function(data, labels){
  pca = c()
  for(label in sort(unique(labels))){
    genenow = which(labels==label)
    if(length(genenow)>1){
      pca_res <- prcomp(t(data[genenow,]), scale. = TRUE)
      pca = rbind(pca, pca_res$x[,1])
    }
    else{
      pca_res <- (data[genenow,]-mean(data[genenow,])) / sd(data[genenow,])
      pca = rbind(pca, pca_res)
    }
  }
  
  #Sigma = cov(t(pca))
  #Omega = 1*(solve(Sigma)!=0)
  Sigma = cor(t(pca))
  rownames(Sigma) = colnames(Sigma) = c(1:nrow(Sigma))
  
  return(list(Sigma = Sigma)) #, Omega=Omega))
}