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


#' Expression Filter 
#' Tests a filter schema on a user defined cut of the data
#' 
#' @usage 
#' 
#' @param exp An expression data set of raw counts
#' @param met A metadata object with columnl lables of; region_label, 
#' external_donor_name_label, broad_cell, cell_type_alias_label, and 
#' sample_name. These correspond to brain region, donor ID, broad cell type 
#' label, specific cell type label, and sample name. 
#' @param pcnt A numeric between 0 and 1 corresponding to a the rank value of the
#' feature's expression to apply the filter of `value`. default = 0.5 (median)
#' @param alue A numeric value to apply the filter to the `pcnt` value
#' @param cell_t character vector of meta data column values to select for broad OR 
#' specific cell type values. If values are broad cell type(s) `is.broad` must
#' be set to TRUE (Optional)
#' @region a character vector of brain region(s)/tissue(s) to subset the data 
#' by (Optional)
#' @param onor_id value of donor ID's to subset data by. (Optional)
#' @param is.broad if cell_t is assigned then is.broad must be set to TRUE 
#' (cell types are broad cell types) or FALSE (cell types are specific cell 
#' types)
#' @param limit Optional. The minimum number of cells allowed to examine for a given cell
#' population. default NULL ie(not set)
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
                        is.broad = NULL, limit = NULL) {
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
    message(paste0("Parameters cannot be combined, one or more of the subsets
                does not exist"))
    return(NA)
  }
  
  # Filter gene expression object by cells
  exp <- as.data.frame(exp[,met$sample_name])
  if(!is.null(limit)){
    if(dim(exp)[1] < limit){
      message(paste0("Cell numbers in ", region, ' - ', cell_t, 
                     ' are below user limit of ', limit ))
      return(NA)
    }
  }
  
  # Filter gene expression object by expression parameters
  initial_counts <- sum(exp)
  
  
  ## Filter out genes of X percent zeros
  present <- function(x, percent = pcnt) mean(x > value) >= percent
  filt_exp <- as.data.frame(exp[apply(exp, 1, present), ])
  
  
  features <- dim(filt_exp)[1]
  count_coverage <- signif( 100*(sum(filt_exp)/initial_counts), 4)
  
  
  return( list(exp = as.data.frame(apply(filt_exp,2,CPM)),
               met = met,
               features = features,
               count_coverage = count_coverage) )
  
}


#' Combined region analysis
#' @param exp gene expression matrix
#' @param met meta data
#' @param region1_name name of region 1 in the combined region analysis
#' @param region2_name name of region 2 in the combined region analysis
#' @return a list object of overall region analysis (list$overall),
#' combined region feature comparison (list$feature_analysis), combined region 
#' count coverage comparison (list$coverage_analysis), and spearman correlation
#' coefficient for the two regions (list$corr)

comb_region <- function(exp, met, region1_name, region2_name){
  
  region1 <- filter_dat(exp = exp, met = met, is.broad = FALSE, 
                              region = region1_name)
  
  region2 <- filter_dat(exp = exp, met = met, is.broad = FALSE, 
                              region = region2_name)
  
  region_comb <- filter_dat(exp = exp, met = met, is.broad = FALSE, 
                            region = c(region1_name, region2_name))
   
  region_overall <- data.frame('subset' = c('combined', 
                                         region1_name, 
                                         region2_name),
                            'features' = c(region_comb$features, 
                                           region1$features, 
                                           region2$features),
                            'count_coverage'= c(region_comb$count_coverage, 
                                                region1$count_coverage, 
                                                region2$count_coverage)
                            )
  
  comb_cell_types <- unique(region_comb$met$cell_type_alias_label)[10:20]
  
  comb_feature_analysis <- region_overall[,1]
  comb_coverage_analysis <- region_overall[,1]
  for (i in comb_cell_types) {
    full <- filter_dat(exp = exp, met = met, region = c(region1_name, region2_name), 
                       is.broad = FALSE, cell_t = i)
    upper <- filter_dat(exp = exp, met = met, region = region1_name, 
                        is.broad = FALSE, cell_t = i)
    lower <- filter_dat(exp = exp, met = met, region = region2_name, 
                        is.broad = FALSE, cell_t = i)
    
    full_f <- full$features
    upper_f <- upper$features
    lower_f <- lower$features
    
    full_c <- full$count_coverage
    upper_c <- upper$count_coverage
    lower_c <- lower$count_coverage
    
    
    comb_feature_analysis <- cbind(comb_feature_analysis, rbind(full_f,upper_f,lower_f))
    comb_coverage_analysis <- cbind(comb_coverage_analysis, rbind(full_c,upper_c,lower_c))
  } 
  
  colnames(comb_feature_analysis) <- c('subset', comb_cell_types)
  colnames(comb_coverage_analysis) <- c('subset', comb_cell_types)
  
  upper_mean_comb <- as.data.frame(apply(upper$exp,1,mean))
  lower_mean_comb <- as.data.frame(apply(lower$exp,1,mean))
  comb_mean <- merge(upper_mean_comb,lower_mean_comb, by =0)
  
  sp_cor <- cor.test(comb_mean[,2], comb_mean[,3], 
                     method = 'spearman', exact = FALSE)$estimate
  #hist(apply(upper$exp,1,median), xlim = c(0,2000), breaks =100)
  #hist(apply(lower$exp,1,median), xlim = c(0,2000), breaks =100)
  
  
  return(list(
    overall = region_overall,
    feature_analysis = comb_feature_analysis,
    coverage_analysis = comb_coverage_analysis,
    corr = sp_cor))
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
  if (result != 'feature_name' && result != 'exp_matrix') {
    stop(paste0('Result must be either feature_name or exp_matrix'))
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
  
    if (is.specific == FALSE) {
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


#' Function to find third quartile CPM value of specific cell types with more
#' than 50 samples by broad cell type
#' 
#'  @param exp gene expression matrix
#'  @param met meta data for expression matrix
#'  @param region region for analysis
#'  @return list object with specific cell type vectors for each broad cell type 


third_quart <- function(exp, met, region) {
  
  ncores <- as.numeric(detectCores() - 1)
  
  # Region subset
  met_region <- met[met$region_label %in% region,]
  
  # Remove any specific cell types with less than 50 cells
  specific_n_50 <- met_region$cell_type_alias_label[table(met_region$cell_type_alias_label) >= 50]
  met_dat <- met_region[met_region$cell_type_alias_label %in% specific_n_50,]
  
  broad_cell_list <- unique(met$broad_cell)[-c(8,9,10)] #not enough Peri, Endo, or VLMC cells
  
  out_list <- c()
  for (i in broad_cell_list) {
    
    # Specific cell list by broad cell type with more than 50 samples
    broad_met <- met_dat[met_dat$broad_cell == i,]
    spec_cell <- unique(broad_met$cell_type_alias_label)
    
    # Specific cell list expression matrix 
    spec_filt <- mclapply(spec_cell, function(x)
      filter_dat(exp = exp, met = broad_met, 
                 pcnt = 0.5, value = 1, 
                 region = region, is.broad = FALSE, 
                 cell_t = x)$exp,
      mc.cores = ncores 
    )
    
    # Function to apply the quantile calculation
    quant_app <- function(x) {return(apply(x,1,function(y) quantile(y,0.75)))}
    
    # Quantile calcs applied to specific cell expression matrices
    spec_quant <- mclapply(spec_filt,mc.cores = ncores, function(x) quant_app(x))
    names(spec_quant) <- spec_cell
    
    command1 <- paste0(
      'out_list <- c(out_list,"', i,'"= list(spec_quant))'
    )
    eval(parse(text=command1))
  }
  
  return(out_list)
}


#' Mean of third quantile list obj function
#' @param input list object of each region i.e. v1c$Unk
#' @return out_mean mean of third quartile for each gene by region

mean_quant_calc <- function(input) {
  out <- c()
  setup <- lapply(input, function(x) rownames_to_column(as.data.frame(x), 
                                                        var = "Row.names"))
  mean_input <- Reduce(function(x,y) merge(x,y, by = 'Row.names', all = TRUE), 
                       setup)
  mean_input <- column_to_rownames(mean_input, var = "Row.names")
  out_mean <- apply(mean_input, 1, function(x) mean(x,na.rm=TRUE))
  return(out_mean)
}


#' List of third quantile mean by region
#' Uses mean_quant_calc to create list objects of cell types by region
#' @param list_obj of third quartile expression for expression by specific
#' cell type, per broad cell type, per region. From third_quant().
#' @return out_list list object of mean third quartile per broad cell type

mean_quant <- function(list_obj) {
  ncores <- as.numeric(detectCores() - 1)
  out_list <-c()
  label <- names(list_obj)
  
  out_list <- c(out_list, mclapply(list_obj, mc.cores = ncores, mean_quant_calc))
  
  names(out_list) <- label
  
  return(out_list)
}


#' Tissue feature analysis
#' 
#' @param list_obj list object of mean expression by features for broad cell type
#' @is.length boolean value indicating whether to return the counts of features or
#' the mean value vector 
#' @return list object for each broad cell type and unique features and features
#' in all cell types

feature_summ <- function(list_obj, is.length = FALSE) {
  
  tot_features <- lapply(list_obj, function(x) unname(names(x)))
  tot_features_unlist <- unlist(tot_features)
  
  overlap <- unname(tot_features_unlist[table(tot_features_unlist) == 9])
  tot_features$tissue_overlap <- unique(overlap)
  
  unique_feat <- unname(tot_features_unlist[table(tot_features_unlist) == 1])
  tot_features$tissue_unique <- unique_feat
  
  if (is.length == TRUE) {
    tot_features <- lapply(tot_features, length)
  }
  
  return(tot_features)
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

