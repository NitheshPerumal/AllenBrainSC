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
  count <- as.data.frame(t(apply(x,1, FUN = median)))
  pruned <- as.data.frame(count[,(count >= z)])
  pruned <- x[names(pruned),]
  return(pruned)
}


#' Summary of Median Cutoffs
#' @param x datframe of gene expression matrix
#' @return summary table of median pruning (using prune function)
#' @return list of total features, and list of features with zeros pruned
#' @return 3 lists of gene features kept with respective cutoffs 
cutoff_summ <- function(x) {
  
  bp <- nrow(x)
  bp_names <<- c(bp_names,rownames(x))
  
  z <- as.data.frame(x[rowSums(x) != 0,])
  z_names <<- c(z_names,rownames(z))
  z <- nrow(z)
  
  cpm_1 <- prune(x,1)
  cpm_1_names <<- c(cpm_1_names, rownames(cpm_1))
  cpm_1 <- nrow(cpm_1)
  
  cpm_0.5 <- prune(x,0.5)
  cpm_0.5_names <<- c(cpm_0.5_names, rownames(cpm_0.5))
  cpm_0.5 <- nrow(cpm_0.5)
  
  cpm_0.1 <- prune(x,0.1)
  cpm_0.1_names <<- c(cpm_0.1_names, rownames(cpm_0.1))
  cpm_0.1 <- nrow(cpm_0.1)
  
  cpm_5 <- prune(x,5)
  cpm_5_names <<- c(cpm_5_names, rownames(cpm_5))
  cpm_5 <- nrow(cpm_5)
  
  cpm_10 <- prune(x,10)
  cpm_10_names <<- c(cpm_10_names, rownames(cpm_10))
  cpm_10 <- nrow(cpm_10)
  
  return(as.data.frame(cbind(bp,z,cpm_1,cpm_0.5,cpm_0.1,cpm_5,cpm_10)))
}


#' Summary of Median Cutoffs for many cells
#' @param x list of objects to pass for cutoff analysis
#' @return summary table of median pruning
cutoff_analysis <- function(cell_list) {
  bp_names <<- c() 
  z_names <<- c()
  cpm_1_names <<- c()
  cpm_0.5_names <<- c()
  cpm_0.1_names <<- c()
  cpm_5_names <<- c()
  cpm_10_names <<- c()
  
  out <- data.frame()
  for (i in cell_list) {
    command3 <- paste0("out <- rbind(out,cutoff_summ(",i,"$exp))")
    eval(parse(text=command3))
  }
  
  unique_features <- cbind(as.numeric(table(table(bp_names))[1]),
                           as.numeric(table(table(z_names))[1]),
                           as.numeric(table(table(cpm_1_names))[1]),
                           as.numeric(table(table(cpm_0.5_names))[1]),
                           as.numeric(table(table(cpm_0.1_names))[1]),
                           as.numeric(table(table(cpm_5_names))[1]),
                           as.numeric(table(table(cpm_10_names))[1]))
  
  cutoff_analysis <- cbind(as.data.frame(cell_list),out)
  cutoff_analysis <- rbind(cutoff_analysis,c('unique_features',unique_features))
  
  rm(bp_names, pos = 1)
  rm(z_names, pos = 1)
  rm(cpm_1_names, pos = 1)
  rm(cpm_0.1_names, pos = 1)
  rm(cpm_0.5_names, pos = 1)
  rm(cpm_10_names, pos = 1)
  rm(cpm_5_names, pos = 1)
  
  return(cutoff_analysis)
}



#' Quantile based pruning
#' Prune features based on quantile cutoff 0.8
#' @param x datframe of gene expression matrix
#' @param z numeric value for threshold of CPM cutoff
#' @return pruned dataframe of gene expression matrix with missing features removed
quant <- function(x,z) {
  count <- as.data.frame(t(as.data.frame(apply(x,1, 
                                               function(x) unname(quantile(x, c(.75))) 
                                               ))))
  pruned <- as.data.frame(count[,count >= z])
  out <- as.data.frame(x[names(pruned),])
  return(out)
}


#' Summary of Quantile Cutoffs for many cells
#' @param x list of objects to pass for quantile analysis
#' @return summary table of quantile pruning
quant_analysis <- function(cell_list) {
  quant_0.1_names <<- c()
  quant_0.5_names <<- c()
  quant_1_names <<- c() 
  quant_5_names <<- c()
  quant_10_names <<- c()
  
  quant_out <- data.frame()
  for (i in cell_list) {
    command3 <- paste0("quant_out <- rbind(quant_out,quant_summ(",i,"$exp))")
    eval(parse(text=command3))
  }
  
  unique_features <- cbind(as.numeric(table(table(quant_1_names))[1]),
                           as.numeric(table(table(quant_0.5_names))[1]),
                           as.numeric(table(table(quant_0.1_names))[1]),
                           as.numeric(table(table(quant_5_names))[1]),
                           as.numeric(table(table(quant_10_names))[1]))
  
  quant_analysis <- cbind(as.data.frame(cell_list), quant_out)
  quant_analysis <- rbind(quant_analysis,c('unique_features', unique_features))
  

  rm(quant_1_names, pos = 1)
  rm(quant_0.1_names, pos = 1)
  rm(quant_0.5_names, pos = 1)
  rm(quant_10_names, pos = 1)
  rm(quant_5_names, pos = 1)
  
  return(quant_analysis)
}


#' Summary of Quantile pruning for different cutoffs
#' @param x datframe of gene expression matrix
#' @return dataframe with summary
#' @return 3 lists of gene features kept with respective cutoffs
quant_summ <- function(x) {
  
  quant_0.1 <- quant(x,0.1)
  quant_0.1_names <<- c(quant_0.1_names, rownames(quant_0.1))
  quant_0.1 <- nrow(quant_0.1)
  
  quant_0.5 <- quant(x,0.5)
  quant_0.5_names <<- c(quant_0.5_names, rownames(quant_0.5))
  quant_0.5 <- nrow(quant_0.5)
  
  quant_1 <- quant(x,1)
  quant_1_names <<- c(quant_1_names, rownames(quant_1))
  quant_1 <- nrow(quant_1)
  
  quant_5 <- quant(x,5)
  quant_5_names <<- c(quant_5_names, rownames(quant_5))
  quant_5 <- nrow(quant_5)
  
  quant_10 <- quant(x,10)
  quant_10_names <<- c(quant_10_names, rownames(quant_10))
  quant_10 <- nrow(quant_10)

  return(as.data.frame(cbind(quant_1,quant_0.5,quant_0.1,quant_5,quant_10) ))
}


#' Function to find summary stats of every feature
#' @param x dataframe of gene matrix to generate summary statistics for
#' @return out dataframe of summary statistics
stats <- function(x) {
  
  wm <- as.data.frame(apply(x, 1, function(x) log2(mean(winsor(x, trim = 0.05, na.rm = TRUE)))))
  colnames(wm)[1] <- 'wins_mean'
  
  m <- as.data.frame(apply(x, 1, function(x) log2(as.numeric(mean(x)))))
  colnames(m)[1] <- 'mean'
  
  med <- as.data.frame(apply(x, 1, function(x) log2(as.numeric(median(x)))))
  colnames(med)[1] <- 'median'
  
  std_dev <- as.data.frame(apply(x, 1, FUN = sd))
  colnames(std_dev)[1] <- 'std_dev'

  feature_name <- as.data.frame(rownames(m))
  colnames(feature_name)[1] <- 'features'
  
  out <- cbind(feature_name, wm, m, med, m-med, std_dev)
  colnames(out)[5] <- 'diff'
  
  hist(out$diff, main = paste0(deparse(substitute(x)),"Mean - Median"), 
       xlab = 'mean - median')
  
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


#' Pruned cell types by region 
#' Prunes features and returns list within list each list is a region with broad
#' broad cell type and specific cell type with feature names and exp matrix
#' 
#' @usage 
#' 
#' @param exp An expression data set of raw counts
#' @param met A metadata object with columnl lables of; region_label, 
#' external_donor_name_label, broad_cell, cell_type_alias_label, and 
#' sample_name. These correspond to brain region, donor ID, broad cell type 
#' label, specific cell type label, and sample name.
#' @region a character vector of brain region(s)/tissue(s) to subset the data 
#' by (Optional)
#' @is.specific boolean to determine if subset with subtypes is wanted, if FALSE
#' broad cell type will be most specific subset, if TRUE specific cell types
#' will be most specific cell type
#' @result c('feature_names','exp_matrix') feature names or expression matrix
#' output determined
#' 
#' @return a list object by region (list$region), broad cell type object 
#' (list$broad cell type), specific cell type object (list$ specific cell type), 
#' feature_name returns list of feature names while exp_matrix returns the 
#' expression matrix for the subset after quantile pruning of 25% 
#' (list$feature_name) or (list$exp_matrix)
#' 
#' @export 
#' 

prune_quant <- function(exp, met, region, is.specific = NULL,
                  result = c('feature_name','exp_matrix')) {
  if (missing(exp)) {
    stop(paste0('Must specify expression matrix'))
  }
  if (missing(met)) {
    stop(paste0('Must specify meta dataframe'))
  }
  if (missing(region)) {
    stop(paste0('Must specify a region'))
  }
  if (missing(is.specific)) {
    stop(paste0('Must specify is.specific with boolean value'))
  }

  broad_cell_list <- unique(met$broad_cell)[-10] #not enough VLMC cells
  
  output_list <- c()
  spec_list <- c()
  if (result == 'feature_name') {
    if(is.specific == FALSE){
      for (i in broad_cell_list) {
        data <- filter_dat(exp = exp, met = met, pcnt = 0.5, value = 1,
                           region = region, is.broad = TRUE, cell_t = i)
        command <- paste0(
          'output_list <- c(output_list,"', i,'"= list(rownames(quant(data$exp,0.5))))'
        )
        eval(parse(text=command))

        }
      }else{

    for (i in broad_cell_list) {
      data <- filter_dat(exp = exp, met = met, pcnt = 0.5, value = 1, 
                         region = region, is.broad = TRUE, cell_t = i)
      
      for (j in unique(data$met$cell_type_alias_label)) {
            
        data_sub <- filter_dat(exp = exp, met = met, pcnt = 0.5, value = 1, 
                              is.broad = FALSE, cell_t = j)
        command <- paste0(
          'spec_list <- c(spec_list,"', j,'"= list(rownames(quant(data_sub$exp,0.5))))'
        )
        eval(parse(text=command))
        }
          
        command1 <- paste0(
        'output_list <- c(output_list,"', i,'"= list(spec_list))'
        )
        eval(parse(text=command1))
      }
    }
  }else{
  
  if(is.specific == FALSE){
    for (i in broad_cell_list) {
      data <- filter_dat(exp = exp, met = met, pcnt = 0.5, value = 1,
                         region = region, is.broad = TRUE, cell_t = i)
      command <- paste0(
        'output_list <- c(output_list,"', i,'"= list(quant(data$exp,0.5)))'
      )
      eval(parse(text=command))
      
    }
  }else{
    
    for (i in broad_cell_list) {
      data <- filter_dat(exp = exp, met = met, pcnt = 0.5, value = 1, 
                         region = region, is.broad = TRUE, cell_t = i)
      
      for (j in unique(data$met$cell_type_alias_label)) {
        
        data_sub <- filter_dat(exp = exp, met = met, pcnt = 0.5, value = 1, 
                               is.broad = FALSE, cell_t = j)
        command <- paste0(
          'spec_list <- c(spec_list,"', j,'"= list(quant(data_sub$exp,0.5)))'
        )
        eval(parse(text=command))
      }
      
      command1 <- paste0(
        'output_list <- c(output_list,"', i,'"= list(spec_list))'
      )
      eval(parse(text=command1))
    }
  }
 }
  
  return(output_list)
  
}


#' Cell Type Specificity Scoring
#' Scores genes based on the specificty to cell types
#' 
#' @usage 
#' 
#' @param cell_list list of object names used for cell type scoring
#' @return a list object of composition table (list$comp) and input data for
#' UpSet plots

composition <- function(cell_list, mainbar_name, sets_name) {
  
  # Dataframe with medians of all features by Broad Cell type
  features <- as.data.frame(rownames(gene_exp))
  colnames(features)[1] <- 'features'
  

  for (i in cell_list) {
    command3 <- paste0("features <- left_join(features,",i,"[, c(1,4)], by = 'features')")
    eval(parse(text=command3))
    #features <- left_join(features, i[, c(1,4)], by = 'features')[,2]
  }
  colnames(features) <- c("features",cell_list)
  
  # Zeroing NA and making column with sum for each gene
  features[is.na(features)] <- 0
  features$sum <- apply(features[,-1], 1, FUN = sum)
  
  
  # Proportion Composition by Median as Cell Type Score
  composition <- features[,2:(ncol(features)-1)]/features$sum
  composition[is.na(composition)] <- 0
  composition <- cbind(features[, 1], composition)
  colnames(composition)[1] <- 'features'
  
  # UpSet plot data
  ups <- features[,2:(ncol(features)-1)]
  ups[ups != 0] <- 1
  
  return( list(comp = composition,
               ups = ups) )
  
}

