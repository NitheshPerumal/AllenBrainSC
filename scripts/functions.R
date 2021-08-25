#' Counts Per Million (CPM) Normalization by sample 
#' @param x datframe of gene expression matrix
#' @return dataframe of normalized expression by Counts per million
CPM <- function(x){
  den <- sum(x)/(10^6)
  x/den
}

#' Median Pruning
#' Removing missing features if 50% or more instances are <= 1 CPM
#' @param x datframe of gene expression matrix
#' @param z numeric value for threshold of CPM cutoff
#' @return pruned dataframe of gene expression matrix with missing features removed
prune <- function(x,z){
  count <- as.data.frame(t(as.data.frame(apply(x[,-2],2, FUN = median))))
  pruned <- as.data.frame(count[,!(count <= z)])
  pruned <- x[,c('sample_name',names(pruned))]
  return(pruned)
}

#' Summary of Median Cutoffs
#' @param x datframe of gene expression matrix
#' @return summary table of median pruning (using prune function) 
cutoff_summ <- function(x){
  y <- x[,-2]
  
  bp <- ncol(y)
  bp_names <<- c(bp_names,names(y))
  
  z <- as.data.frame(y[,colSums(y) != 0])
  z_names <<- c(z_names,names(z))
  z <- ncol(z)
  
  cpm_1 <- prune(x,1)
  cpm_1_names <<- c(cpm_1_names, names(cpm_1))
  cpm_1 <- ncol(cpm_1)
  
  cpm_0.5 <- prune(x,0.5)
  cpm_0.5_names <<- c(cpm_0.5_names, names(cpm_0.5))
  cpm_0.5 <- ncol(cpm_0.5)
  
  cpm_0.1 <- prune(x,0.1)
  cpm_0.1_names <<- c(cpm_0.1_names, names(cpm_0.1))
  cpm_0.1 <- ncol(cpm_0.1)
  
  return(as.data.frame(cbind(bp,z,cpm_1,cpm_0.5,cpm_0.1)))
}

#' Quantile based pruning
#' Prune features based on quantile cutoff 0.8
#' @param x datframe of gene expression matrix
#' @param z numeric value for threshold of CPM cutoff
#' @return pruned dataframe of gene expression matrix with missing features removed
quant <- function(x,z){
  count <- as.data.frame(t(as.data.frame(apply(x[,-2],2, 
                                               function(x) unname(quantile(x, c(.8)))))))
  pruned <- as.data.frame(count[,count >= z])
  out <- as.data.frame(x[,c('sample_name',names(pruned))])
  return(out)
}

#' Summary of Quantile pruning for different cutoffs
#' @param x datframe of gene expression matrix
#' @return dataframe with summary
quant_summ <- function(x){
  
  quant_0.1 <- quant(x,0.1)
  quant_0.1_names <<- c(quant_0.1_names, names(quant_0.1))
  quant_0.1 <- ncol(quant_0.1)
  
  quant_0.5 <- quant(x,0.5)
  quant_0.5_names <<- c(quant_0.5_names, names(quant_0.5))
  quant_0.5 <- ncol(quant_0.5)
  
  quant_1 <- quant(x,1)
  quant_1_names <<- c(quant_1_names, names(quant_1))
  quant_1 <- ncol(quant_1)
  
  
  
  return(as.data.frame(cbind(quant_1,quant_0.5,quant_0.1)))
}

#' Function to find summary stats of every feature
#' @param x dataframe of gene matrix to generate summary statistics for
#' @return out dataframe of summary statistics
stats <- function(x){
  y <- x[,-1]
  
  wm <- as.data.frame(apply(y, 2, function(x) log2(mean(winsor(x, trim = 0.05, na.rm = TRUE)))))
  colnames(wm)[1] <- 'wins_mean'
  
  m <- as.data.frame(apply(y, 2, function(x) log2(as.numeric(mean(x)))))
  colnames(m)[1] <- 'mean'
  
  med <- as.data.frame(apply(y, 2, function(x) log2(as.numeric(median(x)))))
  colnames(med)[1] <- 'median'
  
  std_dev <- as.data.frame(apply(y, 2, FUN = sd))
  colnames(std_dev)[1] <- 'std_dev'
  
  feature_name <- as.data.frame(rownames(m))
  colnames(feature_name)[1] <- 'features'
  
  out <- cbind(feature_name, wm, m, med, m-med, std_dev)
  colnames(out)[5] <- 'diff'
  return(out)
}

