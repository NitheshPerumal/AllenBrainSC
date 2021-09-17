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
m1 <- filter_dat(exp = gene_exp, met = meta, region = c('M1lm','M1ul'))
m1lm <- filter_dat(exp = gene_exp, met = meta, region = 'M1lm')
m1ul <- filter_dat(exp = gene_exp, met = meta, region = 'M1ul')
s1 <- filter_dat(exp = gene_exp, met = meta, region = c('S1lm','S1ul'))
s1ul <- filter_dat(exp = gene_exp, met = meta, region = 'S1ul')
s1lm <- filter_dat(exp = gene_exp, met = meta, region = 'S1lm')
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


# Combined region analysis 
s1_region <- filter_dat(exp = gene_exp, met = meta, is.broad = FALSE, 
                        region = c('S1ul','S1lm'))
s1_comb <- data.frame('subset' = c('s1','s1ul', 's1lm'),
                      'features' = c(s1_region$features, s1ul$features, s1lm$features),
                      'count_coverage'= c(s1_region$count_coverage, s1ul$count_coverage, 
                                          s1lm$count_coverage))

m1_region <- filter_dat(exp = gene_exp, met = meta, is.broad = FALSE, 
                        region = c('M1ul','M1lm'))
m1_comb <- data.frame('subset' = c('m1','m1ul', 'm1lm'),
                      'features' = c(m1_region$features, m1ul$features, m1lm$features),
                      'count_coverage'= c(m1_region$count_coverage, m1ul$count_coverage, 
                                          m1lm$count_coverage))


s1_cell_types <- unique(s1_region$met$cell_type_alias_label)
s1ul_cell_types <- unique(s1ul$met$cell_type_alias_label)
s1lm_cell_types <- unique(s1lm$met$cell_type_alias_label)

s1_subtype <- s1_cell_types[s1_cell_types %in% s1ul_cell_types]
s1_subtype <- s1_subtype[s1_subtype %in% s1lm_cell_types]
s1_subtype <- s1_subtype[c(-5,-23,-31,-33,-34,-43,-73,-78,-85,-90,-99,
                           -101,-102,-104,-105,-107)]

s1_feature_analysis <- s1_comb[1]
s1_coverage_analysis <- s1_comb[1]
for (i in s1_subtype) {
  full <- filter_dat(exp = gene_exp, met = meta, region = c('S1ul','S1lm'), 
                    is.broad = FALSE, cell_t = i)
  upper <- filter_dat(exp = gene_exp, met = meta, region = 'S1ul', 
                     is.broad = FALSE, cell_t = i)
  lower <- filter_dat(exp = gene_exp, met = meta, region = 'S1lm', 
                     is.broad = FALSE, cell_t = i)
  full_f <- full$features
  upper_f <- upper$features
  lower_f <- lower$features
  
  full_c <- full$count_coverage
  upper_c <- upper$count_coverage
  lower_c <- lower$count_coverage
  
  
  s1_feature_analysis <- cbind(s1_feature_analysis, rbind(full_f,upper_f,lower_f))
  s1_coverage_analysis <- cbind(s1_coverage_analysis, rbind(full_c,upper_c,lower_c))
} 
colnames(s1_feature_analysis) <- c('subset', s1_subtype[1:(ncol(s1_feature_analysis)-1)])
colnames(s1_coverage_analysis) <- c('subset', s1_subtype[1:(ncol(s1_feature_analysis)-1)])

upper_mean_s1 <- as.data.frame(apply(upper$exp,1,mean))
lower_mean_s1 <- as.data.frame(apply(lower$exp,1,mean))
s1_mean <- merge(upper_mean_s1,lower_mean_s1, by =0)

cor.test(s1_mean[,2], s1_mean[,3], method = 'spearman', exact = FALSE)
hist(apply(upper$exp,1,median), xlim = c(0,2000), breaks =100, main = 's1 upper')
hist(apply(lower$exp,1,median), xlim = c(0,2000), breaks =100, main = 's1 lower')


# M1 region
m1_cell_types <- unique(m1_region$met$cell_type_alias_label)
m1ul_cell_types <- unique(m1ul$met$cell_type_alias_label)
m1lm_cell_types <- unique(m1lm$met$cell_type_alias_label)

m1_subtype <- m1_cell_types[m1_cell_types %in% m1ul_cell_types]
m1_subtype <- m1_subtype[m1_subtype %in% m1lm_cell_types]
m1_subtype <- m1_subtype[c(-24,-26,-81,-85,-86,-93,-95,-98,-105,-106,
                           -109,-110,-111,-112,-113)]

m1_feature_analysis <- m1_comb[1]
m1_coverage_analysis <- m1_comb[1]
for (i in m1_subtype) {
  full <- filter_dat(exp = gene_exp, met = meta, region = c('M1ul','M1lm'), 
                     is.broad = FALSE, cell_t = i)
  upper <- filter_dat(exp = gene_exp, met = meta, region = 'M1ul', 
                      is.broad = FALSE, cell_t = i)
  lower <- filter_dat(exp = gene_exp, met = meta, region = 'M1lm', 
                      is.broad = FALSE, cell_t = i)
  full_f <- full$features
  upper_f <- upper$features
  lower_f <- lower$features
  
  full_c <- full$count_coverage
  upper_c <- upper$count_coverage
  lower_c <- lower$count_coverage

  
  m1_feature_analysis <- cbind(m1_feature_analysis, rbind(full_f,upper_f,lower_f))
  m1_coverage_analysis <- cbind(m1_coverage_analysis, rbind(full_c,upper_c,lower_c))
} 
colnames(m1_feature_analysis) <- c('subset', m1_subtype[1:(ncol(m1_feature_analysis)-1)])
colnames(m1_coverage_analysis) <- c('subset', m1_subtype[1:(ncol(m1_feature_analysis)-1)])

upper_mean_m1 <- as.data.frame(apply(upper$exp,1,mean))
lower_mean_m1 <- as.data.frame(apply(lower$exp,1,mean))
m1_mean <- merge(upper_mean_m1,lower_mean_m1, by =0)

cor.test(m1_mean[,2], m1_mean[,3], method = 'spearman', exact = FALSE)
hist(apply(upper$exp,1,median), xlim = c(0,2000), breaks =200, main = 'm1 upper')
hist(apply(lower$exp,1,median), xlim = c(0,2000), breaks =200, main = 'm1 lower')


# Heatmap of regions
mtg_mean <- as.data.frame(apply(mtg$exp,1,mean))
a1c_mean <- as.data.frame(apply(a1c$exp,1,mean))
v1c_mean <- as.data.frame(apply(v1c$exp,1,mean))
cgg_mean <- as.data.frame(apply(cgg$exp,1,mean))
s1ul_mean <- as.data.frame(apply(s1ul$exp,1,mean))
s1lm_mean <- as.data.frame(apply(s1lm$exp,1,mean))
m1ul_mean <- as.data.frame(apply(m1ul$exp,1,mean))
m1lm_mean <- as.data.frame(apply(m1lm$exp,1,mean))

region_mean <- transform(merge(mtg_mean,a1c_mean, by = 0, all.x = TRUE), 
                         row.names = Row.names, Row.names = NULL)
region_mean <- transform(merge(region_mean,v1c_mean, by = 0, all.x = TRUE), 
                         row.names = Row.names, Row.names = NULL)
region_mean <- transform(merge(region_mean,cgg_mean, by = 0, all.x = TRUE), 
                         row.names = Row.names, Row.names = NULL)
region_mean <- transform(merge(region_mean,s1ul_mean, by = 0, all.x = TRUE), 
                         row.names = Row.names, Row.names = NULL)
region_mean <- transform(merge(region_mean,s1lm_mean, by = 0, all.x = TRUE), 
                         row.names = Row.names, Row.names = NULL)
region_mean <- transform(merge(region_mean,m1ul_mean, by = 0, all.x = TRUE), 
                         row.names = Row.names, Row.names = NULL)
region_mean <- transform(merge(region_mean,m1lm_mean, by = 0, all.x = TRUE), 
                         row.names = Row.names, Row.names = NULL)
colnames(region_mean) <- c('mtg','a1c','v1c','cgg','s1ul','s1lm','m1ul','m1lm')

region_cor <- cor(region_mean, use = 'complete.obs', method = 'spearman')
heatmap(region_cor, symm = TRUE)

 
# Analyzing different cutoffs to determine cutoff threshold
cell_type_list <- c('exc','inh','unk','astro','oligo','opc','micro')
broad_cell_cutoff_analysis <- cutoff_analysis(cell_type_list)


# Quantile cutoff summary 
cell_type_list <- c('exc','inh','unk','astro','oligo','opc','micro')
broad_cell_quant_analysis <- quant_analysis(cell_type_list)


# List of feature name/exp matrix for cell types by tissue type
mtg_spec_feature_list <- prune_quant(exp = gene_exp, met = meta, region = 'MTG', 
                                     is.specific = TRUE, result ='feature_name')

v1c_spec_feature_list <- prune_quant(exp = gene_exp, met = meta, region = 'V1C', 
                                     is.specific = TRUE, result ='feature_name')

cgg_spec_feature_list <- prune_quant(exp = gene_exp, met = meta, region = 'CgG', 
                                     is.specific = TRUE, result ='feature_name')

m1_spec_feature_list <- prune_quant(exp = gene_exp, met = meta, region = c('M1lm','M1ul'), 
                                     is.specific = TRUE, result ='feature_name')

s1_spec_feature_list <- prune_quant(exp = gene_exp, met = meta, region = c('S1lm','S1ul'), 
                                    is.specific = TRUE, result ='feature_name')

a1c_spec_feature_list <- prune_quant(exp = gene_exp, met = meta, region = 'A1C', 
                                     is.specific = TRUE, result ='feature_name')


# Generating summary statistics 
exc_summary <- stats(exc$exp)
unk_summary <- stats(unk$exp)
inh_summary <- stats(inh$exp)
astro_summary <- stats(astro$exp)
oligo_summary <- stats(oligo$exp)
opc_summary <- stats(opc$exp)
micro_summary <- stats(micro$exp)

mtg_summary <- stats(mtg$exp)
v1c_summary <- stats(v1c$exp)
cgg_summary <- stats(cgg$exp)
m1lm_summary <- stats(m1lm$exp)
s1ul_summary <- stats(s1ul$exp)
s1lm_summary <- stats(s1lm$exp)
m1ul_summary <- stats(m1ul$exp)
a1c_summary <- stats(a1c$exp)


# Cell Type Specificity Scoring ------------------------------------------------

# Broad Cell Type scoring and UpSet plot
broad_cells <- c('exc_summary','inh_summary','unk_summary','astro_summary',
                 'oligo_summary','opc_summary','micro_summary')
broad_cell_comp <- composition(broad_cells)

jpeg(file = '/home/nperumal/AllenBrainSC/plots/UpSet_Broad_Cell_Types.jpeg')
upset(broad_cell_comp$ups, sets = names(broad_cell_comp$ups), 
      order.by = 'freq',
      mainbar.y.label = 'Broad Cell Type Intersection', 
      sets.x.label = 'Broad Cell Type') 
dev.off()


# Scoring by Region and UpSet Plot
region_list <- c('mtg_summary','v1c_summary','cgg_summary','m1lm_summary',
                 's1ul_summary', 's1lm_summary', 'm1ul_summary', 'a1c_summary')
region_comp <- composition(region_list)

jpeg(file = '/home/nperumal/AllenBrainSC/plots/UpSet_Brain_Region.jpeg')
upset(region_comp$ups, 
      sets = names(region_comp$ups), 
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