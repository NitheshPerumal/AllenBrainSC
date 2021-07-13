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
#meta <- read.csv("rstudio/metadata.csv",header = TRUE, sep = ',')
meta <- as.data.frame(data.table::fread(synapser::synGet('syn25883034')$path))

# Gene expression Matrix
#gene_exp <- read.csv("rstudio/matrix.csv", sep = ',')
gene_exp <- as.data.frame(data.table::fread(
    synapser::synGet('syn25883506')$path
  )
)

# Analysis -----------------------------------------------------------------

# Parallel Processing Setup
#ncore <- as.numeric(detectCores()-1)
#cl <- makeCluster(ncore)

# Counts Per Million (CPM) Normalization by sample 
CPM <- function(x){
  den <- sum(x)/(10^6)
  x/den
}

# More efficient but not working
#cpm_exp <- as.data.frame(t(parallel::parApply(cl, gene_exp[,-1], 1, CPM)))
#stopCluster(cl)

# Less Efficient but working
#cpm_exp <- as.data.frame(t(apply(gene_exp[,-1], 1, CPM)))
#cpm_exp <- cbind(gene_exp[,1], cpm_exp)
#colnames(cpm_exp)[1] = "sample_name"

# Getting Normalized data from Synapse
cpm_exp <- as.data.frame(data.table::fread(synapser::synGet('syn25976164')$path))

remove(gene_exp) # Removing un normalized expression matrix to save memory


# Some Exploratory Analysis used to determine how to distinguish cell types
# Counts of Broad Cell Types based on Brain Region
regional_class <- as.data.frame.matrix(xtabs(formula = ~region_label+class_label, meta))
colnames(regional_class)[1] <- "Unlabelled"
write.csv(regional_class, file = 'regional_class_eda.csv', quote = FALSE)

# Counts of Narrow Cell Types based on Brain Region
regional_subclass <- as.data.frame.matrix(xtabs(formula = ~region_label+subclass_label, meta))
colnames(regional_subclass)[1] <- "Unlabelled"
write.csv(regional_subclass, file = 'regional_subclass_eda.csv', quote = FALSE)

# Counts of Narrow Cell Types for Broad Cell Types
class_sub_rel <- as.data.frame.matrix(xtabs(formula = ~subclass_label+class_label, meta))
write.csv(class_sub_rel, file = 'class_sub_rel_eda.csv', quote = FALSE)

# Counts of the most specific distinction vs other two 
cell_type <- as.data.frame.matrix(xtabs(formula = ~cell_type_accession_label+class_label, meta))
cell_type_sub <- as.data.frame.matrix(xtabs(formula = ~cell_type_accession_label+subclass_label, meta))
write.csv(cell_type, file = 'cell_type_eda.csv', quote = FALSE)
write.csv(cell_type_sub, file = 'cell_type_sub_eda.csv', quote = FALSE)


# Broad Cell Type Subsets Metadata
unlabelled <- subset(meta, class_label == "") # Unlabelled samples
inhibitory <- subset(meta, class_label == "GABAergic") # Inhibitory samples
excitatory <- subset(meta, class_label == "Glutamatergic") # Excitatory samples
astro <- subset(meta, subclass_label == 'Astrocyte') # Astrocyte samples
oligo <- subset(meta, subclass_label == 'Oligodendrocyte') # Oligodendrocytes samples
opc <- subset(meta, subclass_label == 'OPC') # Oligo Precursor Cells samples
microglia <- subset(meta, subclass_label == 'Microglia') # Microglia samples
#non_neuronal <- subset(meta, class_label == "Non-neuronal") # Non-neuronal samples


# Subsetting the gene expression matrix by Broad Cell Type
unlab_data <- semi_join(cpm_exp, unlabelled, by = 'sample_name') # Unlabelled Gene exp
inhib_data <- semi_join(cpm_exp, inhibitory, by = 'sample_name') # Inhibitory Gene exp
excit_data <- semi_join(cpm_exp, excitatory, by = 'sample_name') # Excitatory Gene exp
astro_data <- semi_join(cpm_exp, astro, by = 'sample_name')
oligo_data <- semi_join(cpm_exp, oligo, by = 'sample_name')
opc_data <- semi_join(cpm_exp, opc, by = 'sample_name')
microglia_data <- semi_join(cpm_exp, microglia, by = 'sample_name')


# Remvoing un-needed data to free up memory
remove(unlabelled)
remove(inhibitory)
remove(excitatory)
remove(astro)
remove(oligo)
remove(opc)
remove(microglia)


# Removing missing features if 50% or more instances are <= 1 CPM
# @param x datframe of gene expression matrix
# @return pruned dataframe of gene expression matrix with missing features removed
rm_missing <- function(x){
  
  count <- apply(x,2, function(x) sum(x <= 1))
  #count <- parApply(cl,x, 2,function(x) sum(x <= 1))
  pruned <- x[ , -which(names(x) %in% names(which(count >= 0.5*nrow(x))))]
  return(pruned)
}

# Generating dataframe without missing features
#excit_data_rm <- rm_missing(excit_data)
#inhib_data_rm <- rm_missing(inhib_data)
#unlab_data_rm <- rm_missing(unlab_data)
#astro_data_rm <- rm_missing(astro_data)
#oligo_data_rm <- rm_missing(oligo_data)
#opc_data_rm <- rm_missing(opc_data)
#microglia_data_rm <- rm_missing(microglia_data)

# Missing feature pruned data pulling from synapse
excit_data_rm <- as.data.frame(data.table::fread(synapser::synGet('syn25979729')$path))
inhib_data_rm <- as.data.frame(data.table::fread(synapser::synGet('syn25979730')$path))
unlab_data_rm <- as.data.frame(data.table::fread(synapser::synGet('syn25979731')$path))
astro_data_rm <- as.data.frame(data.table::fread(synapser::synGet('syn25979732')$path))
oligo_data_rm <- as.data.frame(data.table::fread(synapser::synGet('syn25979733')$path))
opc_data_rm <- as.data.frame(data.table::fread(synapser::synGet('syn25979734')$path))
microglia_data_rm <- as.data.frame(data.table::fread(synapser::synGet('syn25979735')$path))


# Function to find mean of every column
# @param x dataframe of gene matrix to generate summary statistics for
# @return out dataframe of summary statistics
stats <- function(x){
  y <- x[,c(-1,-2,-3)]
  
  m <- as.data.frame(apply(y, 2, function(x) log2(as.numeric(mean(x)))))
  colnames(m)[1] <- 'mean'
  
  med <- as.data.frame(apply(y, 2, function(x) log2(as.numeric(median(x)))))
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
astro_summary <- stats(astro_data_rm)
oligo_summary <- stats(oligo_data_rm)
opc_summary <- stats(opc_data_rm)
microglia_summary <- stats(microglia_data_rm)


# Histograms of Mean - Median to decide central tendency to use
jpeg(file = '/home/nperumal/AllenBrainSC/plots/Excit_hist.jpeg')
hist(excit_summary$diff, main = 'Excitatory Cells Mean - Median', xlab = 'mean - median')
dev.off()

jpeg(file = '/home/nperumal/AllenBrainSC/plots/Unlab_hist.jpeg')
hist(unlab_summary$diff, main = 'Unlabelled Cells Mean - Median', xlab = 'mean - median')
dev.off()

jpeg(file = '/home/nperumal/AllenBrainSC/plots/Inhib_hist.jpeg')
hist(inhib_summary$diff, main = 'Inhibitory Cells Mean - Median', xlab = 'mean - median')
dev.off()

jpeg(file = '/home/nperumal/AllenBrainSC/plots/Astro_hist.jpeg')
hist(astro_summary$diff, main = 'Astrocytes Mean - Median', xlab = 'mean - median')
dev.off()

jpeg(file = '/home/nperumal/AllenBrainSC/plots/Oligo_hist.jpeg')
hist(oligo_summary$diff, main = 'Oligodendrocytes Mean - Median', xlab = 'mean - median')
dev.off()

jpeg(file = '/home/nperumal/AllenBrainSC/plots/OPC_hist.jpeg')
hist(opc_summary$diff, main = 'Oligo Precrusor Cells Mean - Median', xlab = 'mean - median')
dev.off()

jpeg(file = '/home/nperumal/AllenBrainSC/plots/Microglia_hist.jpeg')
hist(microglia_summary$diff, main = 'Microglia Mean - Median', xlab = 'mean - median')
dev.off()


# Dataframe with medians of all features by Broad Cell type
features <- as.data.frame(names(cpm_exp))
features <- as.data.frame(features[-1,])
colnames(features)[1] <- 'features'

# Median of features by Broad Cell Type
inhib_med <- left_join(features, inhib_summary[,c(1,3)], by = 'features')
unlab_med <- left_join(features, unlab_summary[,c(1,3)], by = 'features')
excit_med <- left_join(features, excit_summary[,c(1,3)], by = 'features')
astro_med <- left_join(features, astro_summary[,c(1,3)], by = 'features')
oligo_med <- left_join(features, oligo_summary[,c(1,3)], by = 'features')
opc_med <- left_join(features, opc_summary[,c(1,3)], by = 'features')
microglia_med <- left_join(features, microglia_summary[,c(1,3)], by = 'features')

features_med <- cbind(inhib_med[,2], unlab_med[,2], excit_med[,2],
                      astro_med[,2], oligo_med[,2], opc_med[,2],
                      microglia_med[,2])
features_med <- cbind(features, features_med)
features_med[is.na(features_med)] <- 0
colnames(features_med) <- c('features', 'Inhibitory','Unlabelled','Excitatory',
                            'Astrocytes', 'Oligodendrocytes', 'OPC', 'Microglia')
features_med$sum <- apply(features_med[,-1],1, FUN = sum)


# Proportion Composition by Median as Cell Type Score
composition <- features_med[,c(-1,-9)]/features_med[,9]
composition[is.na(composition)] <- 0
composition <- cbind(features_med[,1], composition)
colnames(composition)[1] <- 'features'


# Histograms of composition values by cell type
jpeg(file = '/home/nperumal/AllenBrainSC/plots/Inhib_comp.jpeg')
hist(composition$Inhib[composition$Inhib != 0], main = 'Inhibitory Cells Composition values', 
     xlab = 'Composition values')
dev.off()

jpeg(file = '/home/nperumal/AllenBrainSC/plots/Unlab_comp.jpeg')
hist(composition$Unlabelled[composition$Unlabelled != 0], main = 'Unlabelled Cells Composition values', 
     xlab = 'Composition values')
dev.off()

jpeg(file = '/home/nperumal/AllenBrainSC/plots/Excit_comp.jpeg')
hist(composition$Excitatory[composition$Excitatory != 0], main = 'Excitatory Cells Composition values', 
     xlab = 'Composition values')
dev.off()

jpeg(file = '/home/nperumal/AllenBrainSC/plots/Astro_comp.jpeg')
hist(composition$Astrocytes[composition$Astrocytes != 0], main = 'Astrocytes Composition values', 
     xlab = 'Composition values')
dev.off()

jpeg(file = '/home/nperumal/AllenBrainSC/plots/Oligo_comp.jpeg')
hist(composition$Oligodendrocytes[composition$Oligodendrocytes != 0], main = 'Oligodendrocytes Composition values', 
     xlab = 'Composition values')
dev.off()

jpeg(file = '/home/nperumal/AllenBrainSC/plots/OPC_comp.jpeg')
hist(composition$OPC[composition$OPC != 0], main = 'Oligo Precursor Cells Composition values', 
     xlab = 'Composition values')
dev.off()

jpeg(file = '/home/nperumal/AllenBrainSC/plots/Microglia_comp.jpeg')
hist(composition$Microglia[composition$Microglia != 0], main = 'Microglia Composition values', 
     xlab = 'Composition values')
dev.off()


# UpSet Plot
ups <- features_med
ups[ups > 0] <- 1

# Upset plot saved as jpeg
jpeg(file = '/home/nperumal/AllenBrainSC/plots/UpSet_Broad_Cell_Types.jpeg')
upset(ups, sets = c('Inhibitory','Unlabelled','Excitatory', 'Astrocytes',
                    'Oligodendrocytes', 'OPC', 'Microglia'), order.by = 'freq',
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

# Annotations
all.annotations = list(
  dataType = c('clinical','geneExpression'),
  resourceType = 'experimentalData',
  isModelSystem = 'FALSE',
  isMultiSpecimen = 'TRUE',
  fileFormat = 'csv',
  grant = 'na',
  species = 'Human',
  nucleicAcidSource = 'sorted cells',
  organ = 'brain',
  tissue = c('middle temporal gyrus',
             'anterior cingulate cortex',
             'primary visual cortex',
             'primary motor cortex',
             'primary somatosensory cortex',
             'primary auditory cortex'
  ),
  study = c('Allen Institute','Pathway Tracing', 'Treat-AD'), 
  consortium = 'Allen-Brain',
  assay = 'SMART-Seq2'
)

# Counts Table for Broad Cell type per Brain Regions
write.csv(cpm_exp,
          file = 'CPM_Normalized.csv',
          quote = FALSE
)

CPM_dat <- synStore( File(
  path = 'CPM_Normalized.csv',
  name = 'CPM Normalized Gene Expression Matrix',
  parentId = activity$properties$id),
  activityName = activityName,
  activityDescription = activityDescription
)
synapser::synSetAnnotations(CPM_dat, annotations = all.annotations)
file.remove('CPM_Normalized.csv')


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
synapser::synSetAnnotations(CountsTable, annotations = all.annotations)
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
synapser::synSetAnnotations(Inhib, annotations = all.annotations)
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
synapser::synSetAnnotations(Excit, annotations = all.annotations)
file.remove('Excitatory_Cells.csv')


# Non Neuronal Cell Types
#write.csv(non_data,
#          file = 'Non_Neuronal_Cells.csv',
#          quote = FALSE
#)

#Non <- synStore( File(
#  path = 'Non_Neuronal_Cells.csv',
#  name = 'Subset of Non Neuronal Cells from Allen Brain data',
#  parentId = activity$properties$id),
# activityName = activityName,
#  activityDescription = activityDescription
#)
#synapser::synSetAnnotations(Non, annotations = all.annotations)
#file.remove('Non_Neuronal_Cells.csv')


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
synapser::synSetAnnotations(Unlab, annotations = all.annotations)
file.remove('Unlabelled_Cells.csv')


# Missing Features push to synapse (just replaced names for each file)
write.csv(microglia_data_rm,
          file = 'microglia_data_rm.csv',
          quote = FALSE
)

microglia_rm_dat <- synStore( File(
  path = 'microglia_data_rm.csv',
  name = 'Microglia Missing Features pruned',
  parentId = activity$properties$id),
  activityName = activityName,
  activityDescription = activityDescription
)
synapser::synSetAnnotations(microglia_rm_dat, annotations = all.annotations)
file.remove('microglia_data_rm.csv')


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
synapser::synSetAnnotations(Unlab_sum, annotations = all.annotations)
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
synapser::synSetAnnotations(Inhib_sum, annotations = all.annotations)
file.remove('Inhib_Summary.csv')


# Nonneuronal Cells Feature Summary 
#write.csv(non_summary,
#          file = 'Non_Summary.csv',
#         quote = FALSE
#)

#Non_sum <- synStore( File(
#  path = 'Non_Summary.csv',
#  name = 'Non Neuronal Cells Summary Statistics by Feature',
#  parentId = activity$properties$id),
#  activityName = activityName,
#  activityDescription = activityDescription
#)
#synapser::synSetAnnotations(Non_sum, annotations = all.annotations)
#file.remove('Non_Summary.csv')


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
synapser::synSetAnnotations(comp, annotations = all.annotations)
file.remove('composition.csv')


