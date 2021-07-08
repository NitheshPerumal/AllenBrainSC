# Allen Brain Data
# https://portal.brain-map.org/atlases-and-data/rnaseq/human-multiple-cortical-areas-smart-seq


# Libraries ----------------------------------------------------------------

library(tidyverse)
library(synapser)
library(mclust)
library(BiocManager)
library(edgeR)
library(data.table)
library(parallel)
library(UpSetR)

# Load Data -----------------------------------------------------------------

# Metadata for Allen Brain dataset
#--# 1) Added the synapse IDs. It is a pain to download the first time but is
#--# cached after that.
#--# 2) Added in the data.table::fread to quickly load data frames
#--# 3) Formatted the expression to samples=columns and rows=features

#_# meta <- read.csv("rstudio/metadata.csv",header = TRUE, sep = ',')
meta <- data.table::fread(synapser::synGet('syn25883034')$path)

# Gene expression Matrix
#_# gene_exp <- fread("rstudio/matrix.csv", sep = ',')
gene_exp <- t(as.matrix(
  data.table::fread(
    synapser::synGet('syn25883506')$path
  ),
  rownames=1
))


# Analysis -----------------------------------------------------------------

# Parallel Processing Setup
#--# 1) Added the (#') notation to function description and took out log2
#--# 2) Referenced package functions implicity with :: (Will explain in 1 on 1)
#--# 3) Total_cores minus 1 helps to prevent freezes etc.
ncore <- as.numeric(parallel::detectCores()-1)
#--# Make the socket cluster
cl <- parallel::makeCluster(ncore)

# Counts Per Million (CPM) Normalization by sample 
#' @parm x column to perform log2 normalization of CPM
#' @return CPM of variable
CPM <- function(x){
  den <- sum(x)
  x/( den/(10^6) )
}

cpm_exp <- as.data.frame(parallel::parApply(cl, gene_exp, 2, CPM))

# Counts Per Million (CPM) Normalization by sample 
#cpm_exp <- as.data.frame(t(apply(gene_exp[,-1], 1, function(x) (x/sum(x))*1000000)))
#cpm_exp <- cbind(gene_exp[,1], cpm_exp)
#colnames(cpm_exp)[1] = "sample_name"

remove(gene_exp) # Removing un normalized expression matrix to save memory


# Some Exploratory Analysis used to determine how to distinguish cell types
# Counts of Broad Cell Types based on Brain Region
regional_class <- as.data.frame.matrix(xtabs(formula = ~region_label+class_label, meta))
colnames(regional_class)[1] <- "Unlabelled"

# Counts of Narrow Cell Types based on Brain Region
regional_subclass <- as.data.frame.matrix(xtabs(formula = ~region_label+subclass_label, meta))
colnames(regional_subclass)[1] <- "Unlabelled"

# Counts of Narrow Cell Types for Broad Cell Types
class_sub_rel <- as.data.frame.matrix(xtabs(formula = ~subclass_label+class_label, meta))

# Counts of the most specific distinction vs other two 
cell_type <- as.data.frame.matrix(xtabs(formula = ~cell_type_accession_label+class_label, meta))
cell_type_sub <- as.data.frame.matrix(xtabs(formula = ~cell_type_accession_label+subclass_label, meta))


# Broad Cell Type Subsets Metadata
unlabelled <- subset(meta, class_label == "") # Unlabelled samples
inhibitory <- subset(meta, class_label == "GABAergic") # Inhibitory samples
excitatory <- subset(meta, class_label == "Glutamatergic") # Excitatory samples
non_neuronal <- subset(meta, class_label == "Non-neuronal") # Non-neuronal samples

# Double checking that all the data is accounted for
nrow(meta) == nrow(unlabelled) + nrow(inhibitory) + nrow(excitatory) + nrow(non_neuronal)

# Subsetting the gene expression matrix by Broad Cell Type
unlab_data <- semi_join(cpm_exp, unlabelled, by = 'sample_name') # Unlabelled Gene exp
inhib_data <- semi_join(cpm_exp, inhibitory, by = 'sample_name') # Inhibitory Gene exp
excit_data <- semi_join(cpm_exp, excitatory, by = 'sample_name') # Excitatory Gene exp
non_data <- semi_join(cpm_exp, non_neuronal, by = 'sample_name') # Non-neuronal Gene exp

# Remvoing un-needed data to free up memory
remove(unlabelled)
remove(inhibitory)
remove(excitatory)
remove(non_neuronal)


# Removing missing features if 50% or more instances are <= 1 CPM
# @param x datframe of gene expression matrix
# @return pruned dataframe of gene expression matrix with missing features removed
rm_missing <- function(x){
  #count <- apply(x,2, function(x) sum(x <= 1))
  count <- mclapply(x,function(x) sum(x <= 1))
  pruned <- x[ , -which(names(x) %in% names(which(count >= 0.5*nrow(x))))]
  return(pruned)
}


excit_data_rm <- rm_missing(excit_data)
inhib_data_rm <- rm_missing(inhib_data)
unlab_data_rm <- rm_missing(unlab_data)
non_data_rm <- rm_missing(non_data)

# Summary Statistics for each Broad Cell type category
# Function to replace the -Inf with 0
# @param x gene expression dataframe with -Inf values to remove
# @return y dataframe with -Inf replaced with 0
inf <- function(x){
  y <- replace(x, which(x == '-Inf'), as.numeric(0))
  return(y)
}

# Function to find mean of every column
# @param x dataframe of gene matrix to generate summary statistics for
# @return out dataframe of summary statistics
stats <- function(x){
  y <- apply(x[,-1], 2, FUN = inf)
  
  m <- as.data.frame(apply(y, 2, FUN = mean))
  colnames(m)[1] <- 'mean'
  
  med <- as.data.frame(apply(y, 2, FUN = median))
  colnames(med)[1] <- 'median'
  
  std_dev <- as.data.frame(apply(y, 2, FUN = sd))
  colnames(std_dev)[1] <- 'std_dev'
  
  feature_name <- as.data.frame(rownames(m))
  colnames(feature_name)[1] <- 'features'
  
  out <- cbind(feature_name,m,med, m-med, std_dev)
  colnames(out)[4] <- 'diff'
  return(out)
}

excit_summary <- stats(excit_data_rm)
unlab_summary <- stats(unlab_data_rm)
inhib_summary <- stats(inhib_data_rm)
non_summary <- stats(non_data_rm)

hist(excit_summary$diff) # Skewed left
hist(unlab_summary$diff) # Skewed left
hist(inhib_summary$diff) # Skewed left
hist(non_summary$diff) # Looks much more normal

# Dataframe with medians of all features by Broad Cell type
features <- as.data.frame(names(cpm_exp))
features <- as.data.frame(features[-1,])
colnames(features)[1] <- 'features'

# Median of features by Broad Cell Type
inhib_med <- left_join(features, inhib_summary[,c(1,3)], by = 'features')
non_med <- left_join(features, non_summary[,c(1,3)], by = 'features')
unlab_med <- left_join(features, unlab_summary[,c(1,3)], by = 'features')
excit_med <- left_join(features, excit_summary[,c(1,3)], by = 'features')

features_med <- cbind(inhib_med, non_med[,2], unlab_med[,2], excit_med[,2])
features_med[is.na(features_med)] <- 0
colnames(features_med) <- c('features', 'Inhibitory','Nonneuronal','Unlabelled','Excitatory')
features_med$sum <- apply(features_med[,-1],1, FUN = sum)

# Proportion Composition by Median as Cell Type Score
composition <- features_med[,c(-1,-6)]/features_med[,6]
composition[is.na(composition)] <- 0
composition <- cbind(features_med[,1], composition)
colnames(composition)[1] <- 'features'

hist(composition$Inhib[composition$Inhib != 0])
hist(composition$Nonneuronal[composition$Nonneuronal != 0])
hist(composition$Unlabelled[composition$Unlabelled != 0])
hist(composition$Excitatory[composition$Excitatory != 0])

# UpSet Plot
library(UpSetR)
ups <- features_med
ups[ups > 0] <- 1

# Upset plot saved as jpeg
jpeg(file = '/home/nperumal/AllenBrainSC/plots/UpSet_Broad_Cell_Types.jpeg')
upset(ups, sets = c('Inhibitory','Nonneuronal','Unlabelled','Excitatory'), order.by = 'freq',
      mainbar.y.label = 'Broad Cell Type Intersection', sets.x.label = 'Broad Cell Type') 
dev.off()

# Need to make UpSet Plot by region------------------------

# Pushing data to synapse -----------------------------------------------------------

synLogin(authToken = "")

# Setting Synapse ID's
parentID = 'syn25881694'
folder_loc = 'syn25881691'

# Set Activity
activity <- synGet('syn25881691')
activityName = 'Allen Brain Data Analysis'
activityDescription = 'Single Cell analysis of Allen Brain Data'


# Counts Table for Broad Cell type per Brain Regions
write.csv(regional_class,
          file = 'Counts_Broad_Type_vs_Brain_Region.csv',
          quote = FALSE
)

CountsTable <- synStore( File(
  path = 'Counts_Broad_Type_vs_Brain_Region.csv',
  name = 'Counts of Broad Cell Types based on Brain Region',
  parentId = activity$properties$id),
  activityName = activityName,
  activityDescription = activityDescription
)
#synapser::synSetAnnotations(CountsTable, annotations = all.annotations)
file.remove('Counts_Broad_Type_vs_Brain_Region.csv')


# Inhibitory Cell Types
write.csv(inhib_data,
          file = 'Inhibitory_Cells.csv',
          quote = FALSE
)

Inhib <- synStore( File(
  path = 'Inhibitory_Cells.csv',
  name = 'Subset of Inhibitory Cells from Allen Brain data',
  parentId = activity$properties$id),
  activityName = activityName,
  activityDescription = activityDescription
)
#synapser::synSetAnnotations(Inhib, annotations = all.annotations)
file.remove('Inhibitory_Cells.csv')


# Excitatory Cell Types
write.csv(excit_data,
          file = 'Excitatory_Cells.csv',
          quote = FALSE
)

Excit <- synStore( File(
  path = 'Excitatory_Cells.csv',
  name = 'Subset of Excitatory Cells from Allen Brain data',
  parentId = activity$properties$id),
  activityName = activityName,
  activityDescription = activityDescription
)
#synapser::synSetAnnotations(Excit, annotations = all.annotations)
file.remove('Excitatory_Cells.csv')


# Non Neuronal Cell Types
write.csv(non_data,
          file = 'Non_Neuronal_Cells.csv',
          quote = FALSE
)

Non <- synStore( File(
  path = 'Non_Neuronal_Cells.csv',
  name = 'Subset of Non Neuronal Cells from Allen Brain data',
  parentId = activity$properties$id),
  activityName = activityName,
  activityDescription = activityDescription
)
#synapser::synSetAnnotations(Non, annotations = all.annotations)
file.remove('Non_Neuronal_Cells.csv')


# Unlabelled Cells
write.csv(unlab_data,
          file = 'Unlabelled_Cells.csv',
          quote = FALSE
)

Unlab <- synStore( File(
  path = 'Unlabelled_Cells.csv',
  name = 'Subset of Unlabelled Cells from Allen Brain data',
  parentId = activity$properties$id),
  activityName = activityName,
  activityDescription = activityDescription
)
#synapser::synSetAnnotations(Unlab, annotations = all.annotations)
file.remove('Unlabelled_Cells.csv')


# Unlabelled Cells Feature Summary 
write.csv(unlab_summary,
          file = 'Unlabelled_Summary.csv',
          quote = FALSE
)

Unlab_sum <- synStore( File(
  path = 'Unlabelled_Summary.csv',
  name = 'Unlabelled Cells Summary Statistics by Feature',
  parentId = activity$properties$id),
  activityName = activityName,
  activityDescription = activityDescription
)
#synapser::synSetAnnotations(Unlab, annotations = all.annotations)
file.remove('Unlabelled_Summary.csv')


# Inhibitory Cells Feature Summary 
write.csv(inhib_summary,
          file = 'Inhib_Summary.csv',
          quote = FALSE
)

Inhib_sum <- synStore( File(
  path = 'Inhib_Summary.csv',
  name = 'Inhibitory Cells Summary Statistics by Feature',
  parentId = activity$properties$id),
  activityName = activityName,
  activityDescription = activityDescription
)
#synapser::synSetAnnotations(Unlab, annotations = all.annotations)
file.remove('Inhib_Summary.csv')


# Nonneuronal Cells Feature Summary 
write.csv(non_summary,
          file = 'Non_Summary.csv',
          quote = FALSE
)

Non_sum <- synStore( File(
  path = 'Non_Summary.csv',
  name = 'Non Neuronal Cells Summary Statistics by Feature',
  parentId = activity$properties$id),
  activityName = activityName,
  activityDescription = activityDescription
)
#synapser::synSetAnnotations(Unlab, annotations = all.annotations)
file.remove('Non_Summary.csv')


# Percent Composition Data frame
write.csv(composition,
          file = 'composition.csv',
          quote = FALSE
)

comp <- synStore( File(
  path = 'composition.csv',
  name = 'Percent composition of feature by Broad Cell Type',
  parentId = activity$properties$id),
  activityName = activityName,
  activityDescription = activityDescription
)
#synapser::synSetAnnotations(Unlab, annotations = all.annotations)
file.remove('composition.csv')



