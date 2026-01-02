#' Function to correct z-scores by replacing infinite (and large) z-scores by draws from a half-normal distr
#' 
#' @param zscores Vector of z-scores 
#' @param cut Where to cut z-scores, default=NULL in which case set to the largest non-infinite value
#' @param mu Mean of normal to draw from (abs of normal, half-normal), default=NULL in which case set to the largest "non-cut" z-score
#' @param sd SD of half-normal to draw from (abs of normal, half-normal), default=0.5
#' @param pvalues Vector of p-values for ordering, assumed to have same order as zscores
#' @returns Vector of non-infinite corrected z-scores 
rminf <- function(zscores, cut=NULL, mu=NULL, sigma=0.5, pvalues=NULL) {
  
  if (any(is.na(zscores))) {
    stop('NAs in `zscores`, please address and then rerun `rminf`.')
  }
    
  if(sum(zscores == Inf)>0){
    print(paste0('handling ', sum(zscores == Inf), ' infs....'))
    # Sort z-scores to allow for easier replacement, if p-values provided, use those to sort
    sorted_zscores <- sort(zscores, decreasing=TRUE, method='radix')
    if (is.null(pvalues)) {
      sorted_idx_order <- sort(zscores, decreasing=TRUE,index.return=TRUE, method='radix')$ix
    } else {
      if (length(zscores) != length(pvalues)) {
        stop('`zscores` and `pvalues` must have the same length.')
      }
      sorted_idx_order <- sort(pvalues, decreasing=FALSE,index.return=TRUE, method='radix')$ix
    }
    zscores_fixed <- sorted_zscores
    # Define z-scores to be replaced and to keep
    if (is.null(cut)) {
      cut <- max(zscores[is.finite(zscores)])
    }
    zscores_replace <- zscores_fixed[zscores_fixed > cut]
    zscores_keep <- zscores_fixed[zscores_fixed <= cut]
    # If mu is not specified, set to be maximum non-cut z-score
    if (is.null(mu)) {
      mu <- max(zscores_keep)
    }
    # Draw replacement z-scores
    zscores_draws <- mu + abs(rnorm(length(zscores_replace), 0, sigma))
    # Replace z-scores in order (sorted)
    zscores_fixed[zscores_fixed > cut] <- sort(zscores_draws)
    # Reorder fixed z-scores
    zscores_fixed <- zscores_fixed[order(sorted_idx_order)]
    # zscores_fixed <- zscores_fixed[sorted_idx_order]
    if (!is.null(names(zscores))) {
      names(zscores_fixed) <- names(zscores)
    } 
    zscores <- zscores_fixed
  }
  if(sum(zscores == -Inf) > 0){
    print(paste0('handling ', sum(zscores == -Inf), ' -infs....'))
    zscores[which(zscores == -Inf)] <- min(zscores[!is.infinite(zscores)])
  }
  return(zscores)
}

#' Function to plot overlaid histograms
#' 
#' @param vectors_list ...
#' @param breaks ...
#' @param colors ...
#' @param labels ...
#' @param alpha ....
#' @returns Vector of non-infinite corrected z-scores 

overlaid_hists <- function(vectors_list, breaks = 30, colors = NULL, labels = NULL, 
                           alpha = 0.5, main_title = "Overlaid Histograms", 
                           xlab_title = "Value", ylab_title = "Frequency",
                           savepath=NULL, height=4, width=5) {
  # Check if the list is not empty
  if (length(vectors_list) == 0) {
    stop("The list of vectors is empty.")
  }
  
  # If colors are not provided, use a default color palette
  if (is.null(colors)) {
    colors <- RColorBrewer::brewer.pal(max(length(vectors_list), 3), "Set1")
  }
  
  # Create a transparent color vector by using 'rgb' function (alpha channel)
  transparent_colors <- sapply(colors, function(c) rgb(col2rgb(c)[1], col2rgb(c)[2], col2rgb(c)[3], alpha = alpha * 255, maxColorValue = 255))
  
  # Calculate the maximum y-limit across all histograms
  max_y <- 0
  for (i in 1:length(vectors_list)) {
    hist_data <- hist(vectors_list[[i]], breaks = breaks, plot = FALSE)
    max_y <- max(max_y, max(hist_data$counts))
  }
  
  if (!is.null(savepath)) {
    png(savepath, height=height, width=width, units='in', res=300)
  }
  # Create the initial plot for the first histogram
  hist(vectors_list[[1]], col = transparent_colors[1], breaks = breaks, main = main_title, xlab = xlab_title, ylab = ylab_title,
       xlim = c(min(unlist(vectors_list)), max(unlist(vectors_list))), ylim = c(0, max_y), freq = TRUE)
  
  # Add histograms for each additional vector
  for (i in 2:length(vectors_list)) {
    hist(vectors_list[[i]], col = transparent_colors[i], breaks = breaks, add = TRUE, freq = TRUE)
  }
  
  # If labels are not provided, use default labels
  if (is.null(labels)) {
    labels <- paste("Vector", 1:length(vectors_list))
  }
  
  # Add a legend to differentiate the histograms
  legend("topright", legend = labels, fill = transparent_colors[1:length(vectors_list)], cex = 0.8)
  
  if (!is.null(savepath)) {
    dev.off()
  }
}