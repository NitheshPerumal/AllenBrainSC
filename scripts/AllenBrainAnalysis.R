# Allen Brain Data
# https://portal.brain-map.org/atlases-and-data/rnaseq/human-multiple-cortical-areas-smart-seq


# Libraries --------------------------------------------------------------------

library(tidyverse)
library(synapser)
library(mclust)
library(BiocManager)
library(edgeR)
library(data.table)
library(parallel)
library(UpSetR)
library(psych)
source("AllenBrainSC/scripts/functions.R")


# Load Data --------------------------------------------------------------------

# Metadata for Allen Brain dataset
#meta <- read.csv("rstudio/metadata.csv",header = TRUE, sep = ',')
meta <- as.data.frame(
  data.table::fread(
    synapser::synGet('syn25883034')$path
  )
)

# Clean the Metadata
#Remove leading/trailing white space and replace spaces with '_'
meta <- as.data.frame(apply( meta, 2, rep_space ))

#Identify Broad Cell Types
meta$broad_cell <- do.call(
  rbind,
  stringr::str_split(
    meta$cluster_label,
    '_', 
    n = 2, 
    simplify = FALSE
  )
)[,1]

meta$broad_cell[is.na(meta$broad_cell)] <- 'Unk'

meta[ is.na(meta$cell_type_alias_label),]$cell_type_alias_label <- 
  meta[ is.na(meta$cell_type_alias_label),]$outlier_type


# Gene expression Matrix
#gene_exp <- read.csv("rstudio/matrix.csv", sep = ',')
gene_exp <- as.data.frame(
  data.table::fread(
    synapser::synGet('syn25883506')$path
  )
)

# Assign row names
row.names(gene_exp) <- gene_exp$sample_name
gene_exp <- as.data.frame(gene_exp[,2:dim(gene_exp)[2]])
# transpose - assigns gene features to rows, cell IDs to columns
gene_exp <- as.data.frame(t(gene_exp))
# remove gene features with no counts in any cell
gene_exp <- gene_exp[ rowSums(gene_exp != 0) > 0,] 


# Initial EDA ------------------------------------------------------------------

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


# Analysis ---------------------------------------------------------------------

# Operational Definitions used going forward:
# Broad Cell Types (7):
# Excitatory, Inhibitory, Unlabelled, Oligo, OPC, Astro, Microglia

# Specific Cell Types:
# Refer to Broad Cell type + Tissue Location + 2 Markers
# Ex: Astro L1-6 FGFR3 ETNPPL


# Filtering by Broad cell type CPM normalization and Pruning Missing features
exc <- filter_dat(exp = gene_exp, met = meta, is.broad = TRUE, cell_t = 'Exc')
inh <- filter_dat(exp = gene_exp, met = meta, is.broad = TRUE, cell_t = 'Inh')
unk <- filter_dat(exp = gene_exp, met = meta, is.broad = TRUE, cell_t = 'Unk')
astro <- filter_dat(exp = gene_exp, met = meta, is.broad = TRUE, cell_t = 'Astro')
oligo <- filter_dat(exp = gene_exp, met = meta, is.broad = TRUE, cell_t = 'Oligo')
opc <- filter_dat(exp = gene_exp, met = meta, is.broad = TRUE, cell_t = 'OPC')
micro <- filter_dat(exp = gene_exp, met = meta, is.broad = TRUE, cell_t = 'Micro')
endo <- filter_dat(exp = gene_exp, met = meta, is.broad = TRUE, cell_t = 'Endo')
vlmc <- filter_dat(exp = gene_exp, met = meta, is.broad = TRUE, cell_t = 'VLMC')
peri <- filter_dat(exp = gene_exp, met = meta, is.broad = TRUE, cell_t = 'Peri')


# Missing feature pruned data pulling from synapse
# excit_data_rm <- as.data.frame(data.table::fread(
#   synapser::synGet('syn25979729')$path))[, c(-1,-2)]
# inhib_data_rm <- as.data.frame(data.table::fread(
#   synapser::synGet('syn25979730')$path))[, c(-1,-2)]
# unlab_data_rm <- as.data.frame(data.table::fread(
#   synapser::synGet('syn25979731')$path))[, c(-1,-2)]
# astro_data_rm <- as.data.frame(data.table::fread(
#   synapser::synGet('syn25979732')$path))[, c(-1,-2)]
# oligo_data_rm <- as.data.frame(data.table::fread(
#   synapser::synGet('syn25979733')$path))[, c(-1,-2)]
# opc_data_rm <- as.data.frame(data.table::fread(
#   synapser::synGet('syn25979734')$path))[, c(-1,-2)]
# microglia_data_rm <- as.data.frame(data.table::fread(
#   synapser::synGet('syn25979735')$path))[, c(-1,-2)]


# Automatically subset by region and missing features removed
mtg <- filter_dat(exp = gene_exp, met = meta, region = 'MTG')
v1c <- filter_dat(exp = gene_exp, met = meta, region = 'V1C')
cgg <- filter_dat(exp = gene_exp, met = meta, region = 'CgG')
m1lm <- filter_dat(exp = gene_exp, met = meta, region = 'M1lm')
s1ul <- filter_dat(exp = gene_exp, met = meta, region = 'S1ul')
s1lm <- filter_dat(exp = gene_exp, met = meta, region = 'S1lm')
m1ul <- filter_dat(exp = gene_exp, met = meta, region = 'M1ul')
a1c <- filter_dat(exp = gene_exp, met = meta, region = 'A1C')


# Pulling missing feature pruned data for region from synapse
# mtg_data <- as.data.frame(data.table::fread(
#   synapser::synGet('syn25986011')$path))[, -1]
# v1c_data <- as.data.frame(data.table::fread(
#   synapser::synGet('syn25986012')$path))[, c(-1,-2)]
# cgg_data <- as.data.frame(data.table::fread(
#   synapser::synGet('syn25986013')$path))[, c(-1,-2)]
# m1lm_data <- as.data.frame(data.table::fread(
#   synapser::synGet('syn25986014')$path))[, c(-1,-2)]
# s1ul_data <- as.data.frame(data.table::fread(
#   synapser::synGet('syn25986015')$path))[, c(-1,-2)]
# s1lm_data <- as.data.frame(data.table::fread(
#   synapser::synGet('syn25986016')$path))[, c(-1,-2)]
# m1ul_data <- as.data.frame(data.table::fread(
#   synapser::synGet('syn25986017')$path))[, c(-1,-2)]
# a1c_data <- as.data.frame(data.table::fread(
#   synapser::synGet('syn25986018')$path))[, c(-1,-2)]


# Analyzing different cutoffs to determine cutoff threshold
bp_names <<- c() 
z_names <<- c()
cpm_1_names <<- c()
cpm_0.5_names <<- c()
cpm_0.1_names <<- c()

cell_type_list <- c('exc','inh','unk','astro','oligo','opc','micro','endo',
                    'vlmc','peri')

test <- function(x,z) {
  count <- as.data.frame(t(as.data.frame(apply(x[,-2],1, FUN = median))))
  pruned <- as.data.frame(count[!(count <= z)],)
  pruned <- x[,c('sample_name',names(pruned))]
  return(pruned)
}


out <- data.frame()
for (i in cell_type_list) {
  command3 <- paste0("out <- rbind(out,cutoff_summ(",i,"$exp))")
  eval(parse(text=command3))
}

unique_features <- cbind(as.numeric(table(table(bp_names))[1]),
                         as.numeric(table(table(z_names))[1]),
                         as.numeric(table(table(cpm_1_names))[1]),
                         as.numeric(table(table(cpm_0.5_names))[1]),
                         as.numeric(table(table(cpm_0.1_names))[1]))

cutoff_analysis <- cbind(as.data.frame(cell_type_list),out)
cutoff_analysis <- rbind(cutoff_analysis,c('unique_features',unique_features))


# Quantile cutoff summary 
quant_0.1_names <<- c()
quant_0.5_names <<- c()
quant_1_names <<- c()

quant_out <- data.frame()
for (i in cell_type_list) {
  command3 <- paste0("quant_out <- rbind(quant_out,quant_summ(",i,"_data))")
  eval(parse(text=command3))
}

unique_features <- cbind(as.numeric(table(table(quant_1_names))[1]),
                         as.numeric(table(table(quant_0.5_names))[1]),
                         as.numeric(table(table(quant_1_names))[1]))

quant_analysis <- cbind(as.data.frame(cell_type_list),quant_out)
quant_analysis <- rbind(quant_analysis,c('unique_features',unique_features))

# From quant_summ analysis 0.5 has been determined to be good cutoff
# Running quantile pruning and finding number of features with median = 0
zero_list <<- c()
for (i in cell_type_list) {
  command3 <- paste0(i,"_quant<- quant(",i,"_data,0.5)")
  eval(parse(text=command3))
  command5 <- paste0("zero_list<-append(zero_list, 
                     sum(apply(",i,"_quant[,c(-1,-2)],2,FUN=median) == 0))")
  eval(parse(text=command5))
}

# Zero_quant is the summary dataframe with each cell type
# the number of features quant pruned with 0.5 cutoff, the
# number of features with median 0 in quant pruning, and the
# number of features after median pruning with 0.5 cutoff
zero_quant <- cbind(as.data.frame(cell_type_list),as.data.frame(zero_list))
zero_quant <- cbind(zero_quant,quant_analysis$quant_0.5[-8],
                    cutoff_analysis$cpm_0.5[-8])
colnames(zero_quant) <- c('Cell_Type','Median_0','quant_0.5','cpm_0.5')
zero_quant <- zero_quant[, c(1,3,2,4)]



# Nonzero Histograms
for (i in cell_type_list) {
  command4 <- paste0(i,"_hist<-as.data.frame(",i,"_data[,colSums(",i,
                     "_data[, -2]) != 0])")
  eval(parse(text=command4))
}

excit_subs <- apply(excit_hist[, -1], 2, function(x) nonzero_hist(x))
hist(excit_subs, 
     main = 'Excit Nonzero Histogram', xlab = 'log2 percent nonzero')

inhib_subs <- apply(inhib_hist[, -1], 2, function(x) nonzero_hist(x))
hist(inhib_subs, 
     main = 'Inhib Nonzero Histogram', xlab = 'log2 percent nonzero')

unlab_subs <- apply(unlab_hist[, -1], 2, function(x) nonzero_hist(x))
hist(unlab_subs, 
     main = 'Unlabelled Nonzero Histogram', xlab = 'log2 percent nonzero')

astro_subs <- apply(astro_hist[, c(-1,-2)], 2, function(x) nonzero_hist(x))
hist(astro_subs, 
     main = 'Astrocytes Nonzero Histogram', xlab = 'log2 percent nonzero')

oligo_subs <- apply(oligo_hist[, c(-1,-2)], 2, function(x) nonzero_hist(x))
hist(oligo_subs, 
     main = 'Oligodendrocyte Nonzero Histogram', xlab = 'log2 percent nonzero')

opc_subs <- apply(opc_hist[, c(-1,-2)], 2, function(x) nonzero_hist(x))
hist(opc_subs, 
     main = 'OPC Nonzero Histogram', xlab = 'log2 percent nonzero')

microglia_subs <- apply(microglia_hist[, c(-1,-2)], 2, function(x) nonzero_hist(x))
hist(microglia_subs, 
     main = 'Microglia Nonzero Histogram', xlab = 'log2 percent nonzero')


# Generating summary statistics 
excit_summary <- stats(excit_data_rm)
unlab_summary <- stats(unlab_data_rm)
inhib_summary <- stats(inhib_data_rm)
astro_summary <- stats(astro_data_rm)
oligo_summary <- stats(oligo_data_rm)
opc_summary <- stats(opc_data_rm)
microglia_summary <- stats(microglia_data_rm)

mtg_summary <- stats(mtg_data)
v1c_summary <- stats(v1c_data)
cgg_summary <- stats(cgg_data)
m1lm_summary <- stats(m1lm_data)
s1ul_summary <- stats(s1ul_data)
s1lm_summary <- stats(s1lm_data)
m1ul_summary <- stats(m1ul_data)
a1c_summary <- stats(a1c_data)


# Histograms of Mean - Median to decide central tendency to use
jpeg(file = '/home/nperumal/AllenBrainSC/plots/Excit_hist.jpeg')
hist(excit_summary$diff, main = 'Excitatory Cells Mean - Median', 
     xlab = 'mean - median')
dev.off()

jpeg(file = '/home/nperumal/AllenBrainSC/plots/Unlab_hist.jpeg')
hist(unlab_summary$diff, main = 'Unlabelled Cells Mean - Median', 
     xlab = 'mean - median')
dev.off()

jpeg(file = '/home/nperumal/AllenBrainSC/plots/Inhib_hist.jpeg')
hist(inhib_summary$diff, main = 'Inhibitory Cells Mean - Median', 
     xlab = 'mean - median')
dev.off()

jpeg(file = '/home/nperumal/AllenBrainSC/plots/Astro_hist.jpeg')
hist(astro_summary$diff, main = 'Astrocytes Mean - Median', 
     xlab = 'mean - median')
dev.off()

jpeg(file = '/home/nperumal/AllenBrainSC/plots/Oligo_hist.jpeg')
hist(oligo_summary$diff, main = 'Oligodendrocytes Mean - Median', 
     xlab = 'mean - median')
dev.off()

jpeg(file = '/home/nperumal/AllenBrainSC/plots/OPC_hist.jpeg')
hist(opc_summary$diff, main = 'Oligo Precrusor Cells Mean - Median', 
     xlab = 'mean - median')
dev.off()

jpeg(file = '/home/nperumal/AllenBrainSC/plots/Microglia_hist.jpeg')
hist(microglia_summary$diff, main = 'Microglia Mean - Median', 
     xlab = 'mean - median')
dev.off()


# Cell Type Specificity Scoring ------------------------------------------------

# Dataframe with medians of all features by Broad Cell type
features <- as.data.frame(names(cpm_exp))
features <- as.data.frame(features[-1, ])
colnames(features)[1] <- 'features'

# Median of features by Broad Cell Type
inhib_med <- left_join(features, inhib_summary[, c(1,4)], by = 'features')
unlab_med <- left_join(features, unlab_summary[, c(1,4)], by = 'features')
excit_med <- left_join(features, excit_summary[, c(1,4)], by = 'features')
astro_med <- left_join(features, astro_summary[, c(1,4)], by = 'features')
oligo_med <- left_join(features, oligo_summary[, c(1,4)], by = 'features')
opc_med <- left_join(features, opc_summary[, c(1,4)], by = 'features')
microglia_med <- left_join(features, microglia_summary[, c(1,4)], by = 'features')

features_med <- cbind(inhib_med[, 2], unlab_med[, 2], excit_med[, 2],
                      astro_med[, 2], oligo_med[, 2], opc_med[, 2],
                      microglia_med[, 2])
features_med <- cbind(features, features_med)
features_med[is.na(features_med)] <- 0
colnames(features_med) <- c('features', 'Inhibitory','Unlabelled','Excitatory',
                            'Astrocytes', 'Oligodendrocytes', 'OPC', 'Microglia')
features_med[, -1] <- 2^features_med[, -1] # Undo log transform so x > 0
features_med[features_med == 1] <- 0
features_med$sum <- apply(features_med[, -1], 1, FUN = sum)
features_med <- features_med[-1, ]

# Proportion Composition by Median as Cell Type Score
composition <- features_med[, c(-1,-9)]/features_med[, 9]
composition[is.na(composition)] <- 0
composition <- cbind(features_med[, 1], composition)
colnames(composition)[1] <- 'features'


# Histograms of composition values by cell type
jpeg(file = '/home/nperumal/AllenBrainSC/plots/Inhib_comp.jpeg')
hist(composition$Inhib[composition$Inhib != 0], 
     main = 'Inhibitory Cells Composition values', 
     xlab = 'Composition values')
dev.off()

jpeg(file = '/home/nperumal/AllenBrainSC/plots/Unlab_comp.jpeg')
hist(composition$Unlabelled[composition$Unlabelled != 0], 
     main = 'Unlabelled Cells Composition values', 
     xlab = 'Composition values')
dev.off()

jpeg(file = '/home/nperumal/AllenBrainSC/plots/Excit_comp.jpeg')
hist(composition$Excitatory[composition$Excitatory != 0], 
     main = 'Excitatory Cells Composition values', 
     xlab = 'Composition values')
dev.off()

jpeg(file = '/home/nperumal/AllenBrainSC/plots/Astro_comp.jpeg')
hist(composition$Astrocytes[composition$Astrocytes != 0], 
     main = 'Astrocytes Composition values', 
     xlab = 'Composition values')
dev.off()

jpeg(file = '/home/nperumal/AllenBrainSC/plots/Oligo_comp.jpeg')
hist(composition$Oligodendrocytes[composition$Oligodendrocytes != 0], 
     main = 'Oligodendrocytes Composition values', 
     xlab = 'Composition values')
dev.off()

jpeg(file = '/home/nperumal/AllenBrainSC/plots/OPC_comp.jpeg')
hist(composition$OPC[composition$OPC != 0], 
     main = 'Oligo Precursor Cells Composition values', 
     xlab = 'Composition values')
dev.off()

jpeg(file = '/home/nperumal/AllenBrainSC/plots/Microglia_comp.jpeg')
hist(composition$Microglia[composition$Microglia != 0], 
     main = 'Microglia Composition values', 
     xlab = 'Composition values')
dev.off()


# Brain region medians
mtg_med <- left_join(features, mtg_summary[, c(1,4)], by = 'features')
v1c_med <- left_join(features, v1c_summary[, c(1,4)], by = 'features')
cgg_med <- left_join(features, cgg_summary[, c(1,4)], by = 'features')
m1lm_med <- left_join(features, m1lm_summary[, c(1,4)], by = 'features')
s1ul_med <- left_join(features, s1ul_summary[, c(1,4)], by = 'features')
s1lm_med <- left_join(features, s1lm_summary[, c(1,4)], by = 'features')
m1ul_med <- left_join(features, m1ul_summary[, c(1,4)], by = 'features')
a1c_med <- left_join(features, a1c_summary[, c(1,4)], by = 'features')

features_med_region <- cbind(mtg_med[, 2], v1c_med[, 2], cgg_med[, 2],
                      m1lm_med[, 2], s1ul_med[, 2], s1lm_med[, 2],
                      m1ul_med[, 2], a1c_med[, 2])
features_med_region <- cbind(features, features_med_region)
features_med_region[is.na(features_med_region)] <- 0
colnames(features_med_region) <- c('features', 'MTG','V1C','CgG',
                            'M1lm', 'S1ul', 'S1lm', 'M1ul', 'A1C')
features_med_region[, -1] <- 2^features_med_region[, -1] #Undo log transform so x>0
features_med_region[features_med_region == 1] <- 0
features_med_region$sum <- apply(features_med_region[, -1], 1, FUN = sum)
features_med_region <- features_med_region[-1, ]

# Proportion Composition by Median as Cell Type Score
composition_region <- features_med_region[, c(-1,-10)]/features_med_region[, 10]
composition_region[is.na(composition_region)] <- 0
composition_region <- cbind(features_med_region[, 1], composition_region)
colnames(composition_region)[1] <- 'features'


# UpSet Plots
ups <- features_med
ups[ups != 0] <- 1

ups_region <- features_med_region
ups_region[ups_region != 0] <- 1

# Upset plot broad cell type saved as jpeg
jpeg(file = '/home/nperumal/AllenBrainSC/plots/UpSet_Broad_Cell_Types.jpeg')
upset(ups, sets = names(ups)[c(-1,-9)], 
      order.by = 'freq',
      mainbar.y.label = 'Broad Cell Type Intersection', 
      sets.x.label = 'Broad Cell Type') 
dev.off()

#UpSet Plot by region saved as jpeg
jpeg(file = '/home/nperumal/AllenBrainSC/plots/UpSet_Brain_Region.jpeg')
upset(as.data.frame(ups_region), 
      sets = names(ups_region)[c(-1,-10)], 
      order.by = 'freq')
dev.off()


# Pushing data to synapse ------------------------------------------------------

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

