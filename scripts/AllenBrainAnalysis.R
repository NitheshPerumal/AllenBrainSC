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
library(psych)

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


# Removing missing features if 50% or more instances are <= 1 CPM
# @param x datframe of gene expression matrix
# @return pruned dataframe of gene expression matrix with missing features removed
prune <- function(x){
  count <- as.data.frame(t(as.data.frame(apply(x[,-2],2, FUN = median))))
  pruned <- as.data.frame(count[,!(count <= 1)])
  pruned <- x[,c('sample_name',names(pruned))]
  return(pruned)
}


# Missing feature pruning from Broad cell types
broad_type <- as.data.frame(group_by(meta[,c(1,9)], by = 'class_label')[,-3])
for(i in 1:nrow(broad_type)){
  if(broad_type[i,2] == ''){
    broad_type[i,2] <- 'Unlab'
  }
}
for(i in 1:nrow(broad_type)){
  broad_type[i,2] <- switch(broad_type[i,2], 
                          'GABAergic' = 'inhib_data_rm',
                          'Glutamatergic' = 'excit_data_rm',
                          'OPC' = 'opc_data_rm',
                          'Unlab' = 'unlab_data_rm',
                           NA)
}
broad_type <- na.omit(broad_type)
 

for (i in unique(broad_type$class_label)) {
  command <- paste0(i, "<-subset(broad_type, class_label=='", i, "')")
  eval(parse(text=command))
  command2 <- paste0(i, "<-semi_join(cpm_exp,", i,",by = 'sample_name')")
  eval(parse(text=command2))
  command3 <- paste0(i, "<- prune(",i,")")
  eval(parse(text=command3))
}


# Non Neuronal specific cell types feature pruning (decided based on EDA)
sub_type <- as.data.frame(group_by(meta[,c(1,12)], by = 'class_label')[,-3])
for(i in 1:nrow(sub_type)){
  sub_type[i,2] <- switch(sub_type[i,2], 
                      'Astrocyte' = 'astro_data_rm',
                      'Oligodendrocyte' = 'oligo_data_rm',
                      'OPC' = 'opc_data_rm',
                      'Microglia' = 'microglia_data_rm',
                      NA)
}
sub_type <- na.omit(sub_type)
 
 
for (i in unique(sub_type$subclass_label)) {
  command <- paste0(i, "<-subset(sub_type, subclass_label=='", i, "')")
  eval(parse(text=command))
  command2 <- paste0(i, "<-semi_join(cpm_exp,", i,",by = 'sample_name')")
  eval(parse(text=command2))
  command3 <- paste0(i, "<- prune(",i,")")
  eval(parse(text=command3))
} 


# Missing feature pruned data pulling from synapse
excit_data_rm <- as.data.frame(data.table::fread(synapser::synGet('syn25979729')$path))[,c(-1,-2)]
inhib_data_rm <- as.data.frame(data.table::fread(synapser::synGet('syn25979730')$path))[,c(-1,-2)]
unlab_data_rm <- as.data.frame(data.table::fread(synapser::synGet('syn25979731')$path))[,c(-1,-2)]
astro_data_rm <- as.data.frame(data.table::fread(synapser::synGet('syn25979732')$path))[,c(-1,-2)]
oligo_data_rm <- as.data.frame(data.table::fread(synapser::synGet('syn25979733')$path))[,c(-1,-2)]
opc_data_rm <- as.data.frame(data.table::fread(synapser::synGet('syn25979734')$path))[,c(-1,-2)]
microglia_data_rm <- as.data.frame(data.table::fread(synapser::synGet('syn25979735')$path))[,c(-1,-2)]


# Finding total number of unique features
tot_features <- c(names(excit_data_rm), names(inhib_data_rm), names(unlab_data_rm), names(astro_data_rm),
    names(oligo_data_rm), names(opc_data_rm), names(microglia_data_rm))
tot_unique_features <- length(unique(tot_features))


# Automatically subset by region and missing features removed
region <- group_by(meta[,c(1,21)], by = 'region_label')[,-3]
for (i in unique(region$region_label)) {
 command <- paste0(i, "<-subset(region, region_label=='", i, "')")
 eval(parse(text=command))
 command2 <- paste0(i, "<-semi_join(cpm_exp,", i,",by = 'sample_name')")
 eval(parse(text=command2))
 command3 <- paste0(i, "<- prune(",i,")")
 eval(parse(text=command3))
}


# Pulling missing feature pruned data for region from synapse
mtg_rm <- as.data.frame(data.table::fread(synapser::synGet('syn25986011')$path))[,-1]
v1c_rm <- as.data.frame(data.table::fread(synapser::synGet('syn25986012')$path))[,c(-1,-2)]
cgg_rm <- as.data.frame(data.table::fread(synapser::synGet('syn25986013')$path))[,c(-1,-2)]
m1lm_rm <- as.data.frame(data.table::fread(synapser::synGet('syn25986014')$path))[,c(-1,-2)]
s1ul_rm <- as.data.frame(data.table::fread(synapser::synGet('syn25986015')$path))[,c(-1,-2)]
s1lm_rm <- as.data.frame(data.table::fread(synapser::synGet('syn25986016')$path))[,c(-1,-2)]
m1ul_rm <- as.data.frame(data.table::fread(synapser::synGet('syn25986017')$path))[,c(-1,-2)]
a1c_rm <- as.data.frame(data.table::fread(synapser::synGet('syn25986018')$path))[,c(-1,-2)]


# Unlab per region
# unlab <- subset(meta, class_label == '')
# unlab <- unlab[,c(1,21)]
# for(i in 1:nrow(unlab)){
#   unlab[i,2] <- switch(unlab[i,2], 
#                        'MTG' = 'MTG_unlab',
#                        'V1C' = 'V1C_unlab',
#                        'CgG' = 'CgG_unlab',
#                        'M1lm' = 'M1lm_unlab',
#                        'S1ul' = 'S1ul_unlab',
#                        'S1lm' = 'S1lm_unlab',
#                        'M1ul' = 'M1ul_unlab',
#                        'A1C' = 'A1C_unlab')
# }
# 
# for (i in unique(unlab$region_label)) {
#   command <- paste0(i, "<-subset(unlab, region_label=='", i, "')")
#   eval(parse(text=command))
#   command2 <- paste0(i, "<-semi_join(cpm_exp,", i,",by = 'sample_name')")
#   eval(parse(text=command2))
#   command3 <- paste0(i, "<- prune(",i,")")
#   eval(parse(text=command3))
# } 


# Function to find mean of every column
# @param x dataframe of gene matrix to generate summary statistics for
# @return out dataframe of summary statistics
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

excit_summary <- stats(excit_data_rm)
unlab_summary <- stats(unlab_data_rm)
inhib_summary <- stats(inhib_data_rm)
astro_summary <- stats(astro_data_rm)
oligo_summary <- stats(oligo_data_rm)
opc_summary <- stats(opc_data_rm)
microglia_summary <- stats(microglia_data_rm)

# mtg_summary <- stats(mtg_rm)
# v1c_summary <- stats(v1c_rm)
# cgg_summary <- stats(cgg_rm)
# m1lm_summary <- stats(m1lm_rm)
# s1ul_summary <- stats(s1ul_rm)
# s1lm_summary <- stats(s1lm_rm)
# m1ul_summary <- stats(m1ul_rm)
# a1c_summary <- stats(a1c_rm)

mtg_summary <- stats(MTG)
v1c_summary <- stats(V1C)
cgg_summary <- stats(CgG)
m1lm_summary <- stats(M1lm)
s1ul_summary <- stats(S1ul)
s1lm_summary <- stats(S1lm)
m1ul_summary <- stats(M1ul)
a1c_summary <- stats(A1C)


# mtg_unlab_summary <- stats(MTG_unlab)
# v1c_unlab_summary <- stats(V1C_unlab)
# cgg_unlab_summary <- stats(CgG_unlab)
# m1lm_unlab_summary <- stats(M1lm_unlab)
# s1ul_unlab_summary <- stats(S1ul_unlab)
# s1lm_unlab_summary <- stats(S1lm_unlab)
# m1ul_unlab_summary <- stats(M1ul_unlab)
# a1c_unlab_summary <- stats(A1C_unlab)

# length(unique(c(names(MTG_unlab),names(V1C_unlab), names(CgG_unlab), names(M1lm_unlab), names(S1lm_unlab),
#   names(S1ul_unlab), names(M1ul_unlab), names(A1C_unlab))))


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
inhib_med <- left_join(features, inhib_summary[,c(1,4)], by = 'features')
unlab_med <- left_join(features, unlab_summary[,c(1,4)], by = 'features')
excit_med <- left_join(features, excit_summary[,c(1,4)], by = 'features')
astro_med <- left_join(features, astro_summary[,c(1,4)], by = 'features')
oligo_med <- left_join(features, oligo_summary[,c(1,4)], by = 'features')
opc_med <- left_join(features, opc_summary[,c(1,4)], by = 'features')
microglia_med <- left_join(features, microglia_summary[,c(1,4)], by = 'features')

features_med <- cbind(inhib_med[,2], unlab_med[,2], excit_med[,2],
                      astro_med[,2], oligo_med[,2], opc_med[,2],
                      microglia_med[,2])
features_med <- cbind(features, features_med)
features_med[is.na(features_med)] <- 0
colnames(features_med) <- c('features', 'Inhibitory','Unlabelled','Excitatory',
                            'Astrocytes', 'Oligodendrocytes', 'OPC', 'Microglia')
features_med[,-1] <- 2^features_med[,-1] # Undo log transform so x > 0
features_med[features_med == 1] <- 0
features_med$sum <- apply(features_med[,-1], 1, FUN = sum)
features_med <- features_med[-1,]

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


# Brain region medians
mtg_med <- left_join(features, mtg_summary[,c(1,4)], by = 'features')
v1c_med <- left_join(features, v1c_summary[,c(1,4)], by = 'features')
cgg_med <- left_join(features, cgg_summary[,c(1,4)], by = 'features')
m1lm_med <- left_join(features, m1lm_summary[,c(1,4)], by = 'features')
s1ul_med <- left_join(features, s1ul_summary[,c(1,4)], by = 'features')
s1lm_med <- left_join(features, s1lm_summary[,c(1,4)], by = 'features')
m1ul_med <- left_join(features, m1ul_summary[,c(1,4)], by = 'features')
a1c_med <- left_join(features, a1c_summary[,c(1,4)], by = 'features')

features_med_region <- cbind(mtg_med[,2], v1c_med[,2], cgg_med[,2],
                      m1lm_med[,2], s1ul_med[,2], s1lm_med[,2],
                      m1ul_med[,2], a1c_med[,2])
features_med_region <- cbind(features, features_med_region)
features_med_region[is.na(features_med_region)] <- 0
colnames(features_med_region) <- c('features', 'MTG','V1C','CgG',
                            'M1lm', 'S1ul', 'S1lm', 'M1ul', 'A1C')
features_med_region[,-1] <- 2^features_med_region[,-1] # Undo log transform so x > 0
features_med_region[features_med_region == 1] <- 0
features_med_region$sum <- apply(features_med_region[,-1], 1, FUN = sum)
features_med_region <- features_med_region[-1,]

# Proportion Composition by Median as Cell Type Score
composition_region <- features_med_region[,c(-1,-10)]/features_med_region[,10]
composition_region[is.na(composition_region)] <- 0
composition_region <- cbind(features_med_region[,1], composition_region)
colnames(composition_region)[1] <- 'features'


######################################################################
# Brain region medians for Unlabelled subset
# mtg_med <- left_join(features, mtg_unlab_summary[,c(1,4)], by = 'features')
# v1c_med <- left_join(features, v1c_unlab_summary[,c(1,4)], by = 'features')
# cgg_med <- left_join(features, cgg_unlab_summary[,c(1,4)], by = 'features')
# m1lm_med <- left_join(features, m1lm_unlab_summary[,c(1,4)], by = 'features')
# s1ul_med <- left_join(features, s1ul_unlab_summary[,c(1,4)], by = 'features')
# s1lm_med <- left_join(features, s1lm_unlab_summary[,c(1,4)], by = 'features')
# m1ul_med <- left_join(features, m1ul_unlab_summary[,c(1,4)], by = 'features')
# a1c_med <- left_join(features, a1c_unlab_summary[,c(1,4)], by = 'features')
# 
# features_med_region <- cbind(mtg_med[,2], v1c_med[,2], cgg_med[,2],
#                              m1lm_med[,2], s1ul_med[,2], s1lm_med[,2],
#                              m1ul_med[,2], a1c_med[,2])
# features_med_region <- cbind(features, features_med_region)
# features_med_region[is.na(features_med_region)] <- 0
# colnames(features_med_region) <- c('features', 'MTG','V1C','CgG',
#                                    'M1lm', 'S1ul', 'S1lm', 'M1ul', 'A1C')
# features_med_region[,-1] <- 2^features_med_region[,-1] # Undo log transform so x > 0
# features_med_region[features_med_region == 1] <- 0
# features_med_region$sum <- apply(features_med_region[,-1], 1, FUN = sum)
# features_med_region <- features_med_region[-1,]
# 
# # Proportion Composition by Median as Cell Type Score
# composition_region <- features_med_region[,c(-1,-10)]/features_med_region[,10]
# composition_region[is.na(composition_region)] <- 0
# composition_region <- cbind(features_med_region[,1], composition_region)
# colnames(composition_region)[1] <- 'features'
################################################################################


# UpSet Plots
ups <- features_med
ups[ups != 0] <- 1

ups_region <- features_med_region
ups_region[ups_region != 0] <- 1

# Upset plot broad cell type saved as jpeg
jpeg(file = '/home/nperumal/AllenBrainSC/plots/UpSet_Broad_Cell_Types.jpeg')
upset(ups, sets = names(ups)[c(-1,-9)], order.by = 'freq',
      mainbar.y.label = 'Broad Cell Type Intersection', sets.x.label = 'Broad Cell Type') 
dev.off()

#UpSet Plot by region saved as jpeg
jpeg(file = '/home/nperumal/AllenBrainSC/plots/UpSet_Brain_Region.jpeg')
upset(as.data.frame(ups_region), sets = names(ups_region)[c(-1,-10)], order.by = 'freq')
dev.off()


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


# Missing feature brain region push to synapse (just changed the name)
write.csv(MTG,
          file = 'MTG.csv',
          quote = FALSE
)

mtg <- synStore( File(
  path = 'MTG.csv',
  name = 'MTG Missing Features pruned',
  parentId = activity$properties$id),
  activityName = activityName,
  activityDescription = activityDescription
)
synapser::synSetAnnotations(mtg, annotations = all.annotations)
file.remove('MTG.csv')


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

