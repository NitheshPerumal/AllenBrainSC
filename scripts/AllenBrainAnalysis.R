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
cpm_5_names <<- c()
cpm_10_names <<- c()

cell_type_list <- c('exc','inh','unk','astro','oligo','opc','micro')

out <- data.frame()
for (i in cell_type_list) {
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

cutoff_analysis <- cbind(as.data.frame(cell_type_list),out)
cutoff_analysis <- rbind(cutoff_analysis,c('unique_features',unique_features))


# Quantile cutoff summary 
quant_0.1_names <<- c()
quant_0.5_names <<- c()
quant_1_names <<- c()
quant_5_names <<- c()
quant_10_names <<- c()

quant_out <- data.frame()
for (i in cell_type_list) {
  command3 <- paste0("quant_out <- rbind(quant_out,quant_summ(",i,"$exp))")
  eval(parse(text=command3))
}

unique_features <- cbind(as.numeric(table(table(quant_1_names))[1]),
                         as.numeric(table(table(quant_0.5_names))[1]),
                         as.numeric(table(table(quant_1_names))[1]),
                         as.numeric(table(table(quant_5_names))[1]),
                         as.numeric(table(table(quant_10_names))[1]))

quant_analysis <- cbind(as.data.frame(cell_type_list),quant_out)
quant_analysis <- rbind(quant_analysis,c('unique_features',unique_features))


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

# Broad Cell Type scoring
broad_cells <- c('exc_summary','inh_summary','unk_summary','astro_summary',
                 'oligo_summary','opc_summary','micro_summary')
broad_cell_comp <- composition(broad_cells)

jpeg(file = '/home/nperumal/AllenBrainSC/plots/UpSet_Broad_Cell_Types.jpeg')
upset(broad_cell_comp$ups, sets = names(broad_cell_comp$ups)[c(-1,-9)], 
      order.by = 'freq',
      mainbar.y.label = 'Broad Cell Type Intersection', 
      sets.x.label = 'Broad Cell Type') 
dev.off()


# Scoring by Region
region_list <- c('mtg_summary','v1c_summary','cgg_summary','m1lm_summary',
                 's1ul_summary', 's1lm_summary', 'm1ul_summary', 'a1c_summary')
region_comp <- composition(region_list)

jpeg(file = '/home/nperumal/AllenBrainSC/plots/UpSet_Brain_Region.jpeg')
upset(region_comp$ups, 
      sets = names(region_comp$ups)[c(-1,-10)], 
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