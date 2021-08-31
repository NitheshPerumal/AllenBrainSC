#' Counts Per Million (CPM) Normalization by sample
#' @param x datframe of gene expression matrix
#' @return dataframe of normalized expression by Counts per million
CPM <- function(x) {
  den <- sum(x)/(10^6)
  x/den
}


#' Median Pruning
#' Removing missing features if 50% or more instances are <= 1 CPM
#' @param x datframe of gene expression matrix
#' @param z numeric value for threshold of CPM cutoff
#' @return pruned dataframe of gene expression matrix with missing features removed
prune <- function(x,z) {
  count <- as.data.frame(t(as.data.frame(apply(x[,-2],2, FUN = median))))
  pruned <- as.data.frame(count[,!(count <= z)])
  pruned <- x[,c('sample_name',names(pruned))]
  return(pruned)
}


#' Summary of Median Cutoffs
#' @param x datframe of gene expression matrix
#' @return summary table of median pruning (using prune function)
#' @return list of total features, and list of features with zeros pruned
#' @return 3 lists of gene features kept with respective cutoffs 
cutoff_summ <- function(x) {
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
quant <- function(x,z) {
  count <- as.data.frame(t(as.data.frame(apply(x[,-2],2, 
                                               function(x) unname(quantile(x, c(.8)))))))
  pruned <- as.data.frame(count[,count >= z])
  out <- as.data.frame(x[,c('sample_name',names(pruned))])
  return(out)
}


#' Summary of Quantile pruning for different cutoffs
#' @param x datframe of gene expression matrix
#' @return dataframe with summary
#' @return 3 lists of gene features kept with respective cutoffs
quant_summ <- function(x) {
  
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


#' Function to create data for visualization of nonzero histograms
#' @param x dataframe with zero columns removed
#' @return log2 transformation of fraction of nonzero/total samples
nonzero_hist <- function(x) {
  nonzero <- log2(length(x[x != 0])/length(x))
}


#' Function to find summary stats of every feature
#' @param x dataframe of gene matrix to generate summary statistics for
#' @return out dataframe of summary statistics
stats <- function(x) {
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


#' Replace Spaces
#' 
#' This function replaces the spaces in a charcter vector with a user
#' specified character. Specifically ment to be used in an apply loop over
#' a data frame's columns
#' 
#' @export
#' @usage 
#' 
#' @param input a character vector to replace spaces in
#' @param repchar a character to replace the spaces with default = '_'
#' 
#' @examples 
#' example <- c('abc', 'a b c', 'a  b   c', '    ab c ', '', '  ')
#' rep_space( input=example, repchar='_' )
rep_space <- function( input, repchar='_' ){
  #Clean leading spaces
  output <- trimws(input,which = c("both"),whitespace = "[ \t\r\n]")
  
  #Replace in character spaces
  output <- gsub('\\s+',repchar,output)
  
  #Replace empty values with NA
  output <- dplyr::na_if(output, '')
  
  return(output)
}


#' Expression Filter Test
#' Tests a filter schema on a user defined cut of the data
#' 
#' @usage 
#' 
#' @param exp An expression data set of raw counts
#' @param met A metadata object with columnl lables of; region_label, 
#' external_donor_name_label, broad_cell, cell_type_alias_label, and 
#' sample_name. These correspond to brain region, donor ID, broad cell type 
#' label, specific cell type label, and sample name. 
#' @pcnt A numeric between 0 and 1 corresponding to a the rank value of the
#' feature's expression to apply the filter of `value`. default = 0.5 (median)
#' @value A numeric value to apply the filter to the `pcnt` value
#' @cell_t character vector of meta data column values to select for broad OR 
#' specific cell type values. If values are broad cell type(s) `is.broad` must
#' be set to TRUE (Optional)
#' @region a character vector of brain region(s)/tissue(s) to subset the data 
#' by (Optional)
#' @donor_id value of donor ID's to subset data by. (Optional)
#' @is.broad if cell_t is assigned then is.broad must be set to TRUE 
#' (cell types are broad cell types) or FALSE (cell types are specific cell 
#' types)
#' 
#' @return a list object of filtered counts matrix with CPM normalization 
#' (list$exp), filtered meta data (list$met), number of features remaining 
#' (list$features), and percent of counts retained after filtering counts 
#' based on percentile cutoff in the filtered cell population (list$count_coverage)
#' 
#' @export 
#' 
#' @examples 
#'set.seed(42)
#'s_meta <- as.data.frame(cbind(
#'  sample_name = sample(random::randomStrings(
#'    n=1000, len=6, digits=F, upperalpha=TRUE, loweralpha=F, unique=T), 1000, 
#'    replace = F),
#'  cell_type_alias_label = sample( random::randomStrings(
#'    n=17, len=6, digits=F, upperalpha=TRUE, loweralpha=F, unique=F), 1000, 
#'    replace = T),
#'  region_label = sample( random::randomStrings(
#'    n=7, len=3, digits=F, upperalpha=TRUE, loweralpha=F, unique=F), 1000, 
#'    replace = T),
#'  external_donor_name_label = sample( random::randomStrings(
#'    n=4, len=6, digits=F, upperalpha=TRUE, loweralpha=F, unique=F), 1000, 
#'    replace = T),
#'  broad_cell = sample( random::randomStrings(
#'    n=7, len=5, digits=F, upperalpha=TRUE, loweralpha=F, unique=F), 1000, 
#'    replace = T)
#'))
#'
#'s_exp <- (rbind(
#'  rnbinom(n=1000, 20, p=.75 ),
#'  rnbinom(n=1000, 20, p=.95 ),
#'  rnbinom(n=1000, 20, p=.25 ),
#'  rnbinom(n=1000, 4, p=.25 ),
#'  rnbinom(n=1000, 4, p=.95 ),
#'  rnbinom(n=1000, 4, p=.75 ),
#'  rnbinom(n=1000, 8, p=.25 ),
#'  rnbinom(n=1000, 8, p=.95 ),
#'  rnbinom(n=1000, 8, p=.75 ),
#'  rnbinom(n=1000, 14, p=.5 )
#'))
#'s_exp <- s_exp[rep(1:nrow(s_exp), times = 100), ]
#'colnames(s_exp) <- s_meta$sample_name
#'example <- filter_test(met = s_meta, exp = s_exp, pcnt = .9, value = 19)

filter_dat <- function(exp, met, pcnt = .5, value = 1, 
                        cell_t = NULL, donor = NULL, region = NULL, 
                        is.broad = NULL) {
  if(missing(cell_t)){
    cell_t <- NULL
  }
  if(missing(donor)){
    donor <- NULL
  }
  if(missing(region)){
    region <- NULL
  }
  if(missing(is.broad)){
    is.broad <- NULL
  }
  #Confirm is.broad is set correctly if cell_t is specified
  if(is.null(cell_t)) {
  }else{
    if(!is.null(cell_t) & is.null(is.broad)) {
      stop(
        paste0('Must specify is.broad as TRUE or FALSE if cell_t is specified')
      )
    }else{
      if(!is.null(cell_t) & (isTRUE(is.broad) | isFALSE(is.broad))) {
        
      }else{
        stop(
          paste0('is.broad must be set to TRUE or FALSE')
        )
      }
    }
  }
  # Filter meta data
  if(!is.null(region)){
    met <- met[ met$region_label %in% region, ]
  }
  if(!is.null(donor)){
    met <- met[ met$external_donor_name_label %in% donor, ]
  }
  if(!is.null(cell_t)) {
    if(is.broad == TRUE){
      met <- met[ met$broad_cell %in% cell_t, ]
    }else{
      met <- met[ met$cell_type_alias_label %in% cell_t, ]
    }
  }

  # In case the combination of parameters cannot be used to subset
  if(nrow(met) == 0){
    stop(paste0("Parameters cannot be combined, one or more of the subsets
                does not exist"))
  }
  
  # Filter gene expression object by cells
  exp <- as.data.frame(exp[,met$sample_name])
  
  
  # Filter gene expression object by expression parameters
  initial_counts <- sum(exp)
  
  
  ## Filter out genes of X percent zeros
  present <- function(x, percent = pcnt) mean(x > value) >= percent
  filt_exp <- exp[apply(exp, 1, present), ]
  
  
  features <- dim(filt_exp)[1]
  count_coverage <- signif( 100*(sum(filt_exp)/initial_counts), 4)
  
  
  return( list(exp = as.data.frame(apply(filt_exp,2,CPM)), 
               met = met, 
               features = features, 
               count_coverage = count_coverage) )
  
}

