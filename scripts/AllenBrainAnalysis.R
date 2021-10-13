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

source("scripts/functions.R")

 
# Load Data --------------------------------------------------------------------

#_# Add Login Step
synapser::synLogin()

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

#_# cell_type_alias_label is the more specific cell type
# Counts of Narrow Cell Types for Broad Cell Types
class_sub_rel_2 <- as.data.frame.matrix(xtabs(formula = ~cell_type_alias_label+class_label, meta))
write.csv(class_sub_rel_2, file = 'class_sub_rel_2_eda.csv', quote = FALSE)

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

#_# Move combination Up Here:
#Combined region analysis for s1 and m1
s1_comb_analysis <- comb_region(exp = gene_exp, met = meta, 'S1ul', 'S1lm')
m1_comb_analysis <- comb_region(exp = gene_exp, met = meta, 'M1ul', 'M1lm')

table(meta$region_label)
meta[ grepl('S1',meta$region_label), ]$region_label <- 'S1'
meta[ grepl('M1',meta$region_label), ]$region_label <- 'M1'
table(meta$region_label)

################################################################################
#_# I'm going to test the per-region filtering compared to not.
#_# I'm going to use micro-glia becuase there's no sub-type substructure
micro <- filter_dat(exp = gene_exp, met = meta, is.broad = TRUE, cell_t = 'Micro')
micro_sink <- micro
dim(micro_sink$exp)
# 796 750
micro <- filter_dat(exp = gene_exp, met = meta, is.broad = TRUE, pcnt = .25, value = .1, cell_t = 'Micro')
dim(micro$exp)
# 5436 750


micro_test <- list(A1C=NULL,
                   CgG=NULL,
                   M1=NULL,
                   MTG=NULL,
                   S1=NULL,
                   V1C=NULL
)
features <- list()
total_features <- NULL
for( region in names(micro_test)) {
  micro_test[[region]] <- filter_dat(
    exp = gene_exp, met = meta,
    is.broad = FALSE,
    pcnt = .25, value = .1, 
    cell_t = 'Micro_L1-6_C1QC', region = region
)
  features[[region]] <- row.names(micro_test[[region]]$exp)
  total_features <- c(total_features,row.names(micro_test[[region]]$exp))
  total_features <- total_features[!duplicated(total_features)]
}
ComplexHeatmap::UpSet( ComplexHeatmap::make_comb_mat(features),
                       comb_order = order(ComplexHeatmap::comb_size(
                         ComplexHeatmap::make_comb_mat(features)
                       )
                      )
)
length(total_features)
# 7408 unique features across 6 tissues
dim(micro$exp)
# 5436 Features if we don't look by region

#_# Comp to sorted cells from validation data:
vali <- read.table(synapser::synGet('syn26275022')$path, 
                          header = T, sep = '\t'
                  ) %>%
  tibble::column_to_rownames('feature') 
vali <- as.data.frame(edgeR::cpm(as.matrix(vali)))

cov <- read.table(synapser::synGet('syn26275021')$path, 
                  header = T, sep = '\t'
) %>%
  tibble::column_to_rownames('sample')

vali_filt <- vali[ , row.names(cov[cov$cell_type == 'Myeloid',])]
vali_filt <- apply(vali_filt, 1, mean)

table(vali_filt>=1) 
# FALSE  TRUE 
# 11030 14343

#load biomart
bm <- read.table(synapser::synGet('syn26275023')$path, 
                   header = T, sep = '\t'
) %>%
  tibble::column_to_rownames('ensembl_gene_id')

table( total_features %in% bm[names(vali_filt>=1),]$hgnc_symbol )
#FALSE  TRUE 
#1202  6206
table( row.names(micro$exp) %in% bm[names(vali_filt>=1),]$hgnc_symbol )
################################################################################

######################
#_# look into first filtering along fine cell type with-in region THEN combining
Exc <- list(A1C=list(),
            CgG=list(),
            M1=list(),
            MTG=list(),
            S1=list(),
            V1C=list()
)
features_exc <- list()
total_features_exc <- NULL

cells <- names(table(meta[meta$broad_cell == 'Exc',]$cell_type_alias_label))
# Loop Regions names(Exc)[]
for(region in names(Exc)) {
  # Loop Specific Cell_types
  for(cell_type in cells) {
    Exc[[region]][[cell_type]] <- filter_dat(
      exp = gene_exp, met = meta,
      is.broad = FALSE,
      pcnt = .5, value = 1, 
      cell_t = cell_type, region = region, limit = 50
    )
  }
}

# UpSet Plot of overlap by Tissue Type
Exc_tissues <- list(A1C=NULL,
                    CgG=NULL,
                    M1=NULL,
                    MTG=NULL,
                    S1=NULL,
                    V1C=NULL
)
for(region in names(Exc)) {
  # Loop Specific Cell_types
  for(cell_type in cells) {
    if( !is.na(Exc[[region]][[cell_type]]) ){
      Exc_tissues[[region]] <- c(Exc_tissues[[region]], 
                                 row.names(Exc[[region]][[cell_type]]$exp)
      )
      Exc_tissues[[region]] <- 
        Exc_tissues[[region]][!duplicated(Exc_tissues[[region]] )]
    }
  }
}
ComplexHeatmap::UpSet( ComplexHeatmap::make_comb_mat(Exc_tissues),
                       comb_order = order(ComplexHeatmap::comb_size(
                         ComplexHeatmap::make_comb_mat(Exc_tissues)
                       )
                       )
)

# UpSet Plot of overlap by Cell Types
for(celltype in cells) {
  Exc_cell_types <- list(celltype=NULL)
}

for(celltype in cells) {
  # Loop Specific Cell_types
  for(region in names(Exc)) {
    if( !is.na(Exc[[region]][[celltype]]) ){
      Exc_cell_types[[celltype]] <- c(Exc_cell_types[[celltype]], 
                                 row.names(Exc[[region]][[celltype]]$exp)
      )
      Exc_cell_types[[celltype]] <- 
        Exc_cell_types[[celltype]][!duplicated(Exc_cell_types[[celltype]] )]
    }
  }
}

### Make a heatmap of overlap
genedf <- stack(setNames(Exc_cell_types, nm=names(Exc_cell_types)))
heatmap_df <- table(genedf[2:1]) %*% t(table(genedf[2:1]))
heatmap_df <- as.data.frame(heatmap_df)

## normalize to percent overlap
for(row in 1:dim(heatmap_df)[1]){
  for(col in 1:dim(heatmap_df)[1]){
    if(row == col){
      #if on diagonal value will be one anyways, take row
      heatmap_df[row,col] <- heatmap_df[row,col] / 
        length(Exc_cell_types[[row.names(heatmap_df)[row]]])
    }else{
      if( row < col ){
        #closer to row name then normalize by that value
        heatmap_df[row,col] <- heatmap_df[row,col] / 
          length(Exc_cell_types[[row.names(heatmap_df)[row]]])
      }else{
        #closer to col name then normalize by that value
        heatmap_df[row,col] <- heatmap_df[row,col] / 
          length(Exc_cell_types[[colnames(heatmap_df)[col]]])
      }
    }
  }
}

ComplexHeatmap::Heatmap(heatmap_df)

## Compare to the validation dataset
vali_filt <- row.names(vali[ vali[ , row.names(cov[cov$cell_type == 'Neuron',])] >1,
])

table( bm[vali_filt,]$hgnc_symbol %in% Exc_tissues$A1C )
table( bm[vali_filt,]$hgnc_symbol %in% Exc_tissues$CgG )
table( bm[vali_filt,]$hgnc_symbol %in% Exc_tissues$M1 )
table( bm[vali_filt,]$hgnc_symbol %in% Exc_tissues$MTG )
table( bm[vali_filt,]$hgnc_symbol %in% Exc_tissues$S1 )
table( bm[vali_filt,]$hgnc_symbol %in% Exc_tissues$V1C )

total_features_exc <- c(Exc_tissues$A1C, Exc_tissues$CgG, Exc_tissues$M1,
                        Exc_tissues$MTG, Exc_tissues$S1, Exc_tissues$V1C, 
                        Exc_tissues$CgG
                      )
total_features_exc <- total_features_exc[!duplicated(total_features_exc)]
table( bm[vali_filt,]$hgnc_symbol %in% total_features_exc )


#FALSE  TRUE 
#1202  6206
table( row.names(micro$exp) %in% bm[names(vali_filt>=1),]$hgnc_symbol )

######################

# Filtering by Broad cell type CPM normalization and Pruning Missing features
exc <- filter_dat(exp = gene_exp, met = meta, is.broad = TRUE, cell_t = 'Exc')
inh <- filter_dat(exp = gene_exp, met = meta, is.broad = TRUE, cell_t = 'Inh')
unk <- filter_dat(exp = gene_exp, met = meta, is.broad = TRUE, cell_t = 'Unk')
astro <- filter_dat(exp = gene_exp, met = meta, is.broad = TRUE, cell_t = 'Astro')
oligo <- filter_dat(exp = gene_exp, met = meta, is.broad = TRUE, cell_t = 'Oligo')
opc <- filter_dat(exp = gene_exp, met = meta, is.broad = TRUE, cell_t = 'OPC')
micro <- filter_dat(exp = gene_exp, met = meta, is.broad = TRUE, cell_t = 'Micro')
endo <- filter_dat(exp = gene_exp, met = meta, is.broad = TRUE, cell_t = 'Endo')
peri <- filter_dat(exp = gene_exp, met = meta, is.broad = TRUE, cell_t = 'Peri')

#_# This can be a named list obj
#broad_feature_list <-c(rownames(exc$exp), rownames(inh$exp), rownames(unk$exp),
#                       rownames(astro$exp), rownames(oligo$exp), rownames(opc$exp),
#                       rownames(micro$exp), rownames(endo$exp), rownames(peri$exp))
broad_feature_list <- list(
  'Exc' = rownames(exc$exp), "Inh" = rownames(inh$exp), 'Unk' = rownames(unk$exp), 
  'Astro' = rownames(astro$exp), 'Oligo' = rownames(oligo$exp), 
  'OPC' = rownames(opc$exp),'Micro' = rownames(micro$exp), 
  'Endo' = rownames(endo$exp),'Peri' = rownames(peri$exp)
)
#broad_feature_unique <- table(table(broad_feature_list))[1]
broad_feature_unique <- as.numeric(table(table(unlist(broad_feature_list)))[1])
#broad_feature_overlap <- table(table(broad_feature_list))[9]
broad_feature_overlap <- as.numeric(table(table(unlist(broad_feature_list)))[9])

#broad_feature_count <- as.data.frame(c(exc$features, inh$features, unk$features, 
#                                       astro$features, oligo$features, 
#                                       opc$features, micro$features, 
#                                       endo$features, peri$features,
#                                       length(unique(broad_feature_list)) ))

#broad_cell_list <- as.data.frame(c(unique(meta$broad_cell)[-10], 'all_cells'))
#feature_overview <- cbind(broad_cell_list, broad_feature_count)
#colnames(feature_overview) <- c('Cell type','Feature count')

#_# Try not to use spaces it will really mess with code later on
feature_overview <- data.frame(`Cell type` = names(broad_feature_count),
                                  `Feature count` = as.character(unlist(broad_feature_count))
                                )

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
#_# Given the amount of heterogenity across cell type we should re-visit this another way in a while
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

 
# Analyzing different Median cutoffs to determine cutoff threshold
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


# Third quartile and mean of third quart for each broad cell type
mtg_third_quart <- third_quart(exp = gene_exp, met = meta, region = 'MTG') 
v1c_third_quart <- third_quart(exp = gene_exp, met = meta, region = 'V1C')
cgg_third_quart <- third_quart(exp = gene_exp, met = meta, region = 'CgG')
m1_third_quart <- third_quart(exp = gene_exp, met = meta, region = c('M1lm','M1ul'))
s1_third_quart <- third_quart(exp = gene_exp, met = meta, region = c('S1lm','S1ul'))
a1c_third_quart <- third_quart(exp = gene_exp, met = meta, region = 'A1C')


mtg_mean_quant <- mean_quant(mtg_third_quart)
v1c_mean_quant <- mean_quant(v1c_third_quart)
cgg_mean_quant <- mean_quant(cgg_third_quart)
m1_mean_quant <- mean_quant(m1_third_quart)
s1_mean_quant <- mean_quant(s1_third_quart)
a1c_mean_quant <- mean_quant(a1c_third_quart)

# mtg_feature_summ <- feature_summ(mtg_mean_quant, is.length = TRUE)
# v1c_feature_summ <- feature_summ(v1c_mean_quant, is.length = TRUE)
# cgg_feature_summ <- feature_summ(cgg_mean_quant, is.length = TRUE)
# m1_feature_summ <- feature_summ(m1_mean_quant, is.length = TRUE)
# s1_feature_summ <- feature_summ(s1_mean_quant, is.length = TRUE)
# a1c_feature_summ <- feature_summ(a1c_mean_quant, is.length = TRUE)
# 
# tissue_feature_summ <- rbind(as.data.frame(mtg_feature_summ), as.data.frame(v1c_feature_summ),
#                              as.data.frame(cgg_feature_summ), as.data.frame(m1_feature_summ),
#                              as.data.frame(s1_feature_summ), as.data.frame(a1c_feature_summ))
# 
# rownames(tissue_feature_summ) <- c('mtg','v1c','cgg','m1','s1','a1c')


mtg_feature_summ <- feature_summ(mtg_mean_quant)
v1c_feature_summ <- feature_summ(v1c_mean_quant)
cgg_feature_summ <- feature_summ(cgg_mean_quant)
m1_feature_summ <- feature_summ(m1_mean_quant)
s1_feature_summ <- feature_summ(s1_mean_quant)
a1c_feature_summ <- feature_summ(a1c_mean_quant)


ups_feature_summ <- as.data.frame(rownames(gene_exp))
ups_feature_summ <- column_to_rownames(ups_feature_summ, var = 'rownames(gene_exp)')
for (i in unique(meta$broad_cell)[-10]) {
  command <- paste0('hist_list <- c(mtg_feature_summ$',i, ',v1c_feature_summ$',i,
                    ',cgg_feature_summ$',i, ',m1_feature_summ$',i,
                    ',s1_feature_summ$',i, ',a1c_feature_summ$',i,')')
  eval(parse(text=command))
  hist_in <- table(table(hist_list))
  
  hist_frame <- as.data.frame(unique(hist_list))
  hist_frame <- column_to_rownames(hist_frame, var = 'unique(hist_list)')
  hist_frame[,1] <- 1
  ups_feature_summ <- merge(ups_feature_summ, hist_frame, by = 0, all = TRUE)
  ups_feature_summ <- column_to_rownames(ups_feature_summ, var = 'Row.names')
  colnames(ups_feature_summ)[ncol(ups_feature_summ)] <- i
  
  h <- barplot(hist_in, main = i)
  text(h, 4000 , paste(hist_in, sep="") ,cex=1)
}

ups_feature_summ[is.na(ups_feature_summ)] <- 0
upset(ups_feature_summ, sets = names(ups_feature_summ), 
      order.by = 'freq',
      mainbar.y.label = 'Feature Overlap by Cells', 
      sets.x.label = 'Broad Cell Type') 

ups_feature_summ <- as.data.frame(ups_feature_summ)

#ups_13_feat <- rownames(ups_feature_summ[apply(ups_feature_summ, 1, sum) == 4,])


cells <- c('Endo', 'Micro', 'Peri', 'OPC')
for (i in cells) {
    
  command <- paste0('
  hist_dat <- cbind(
    as.data.frame(mtg_mean_quant$',i,'[ups_13_feat], ups_13_feat),
    as.data.frame(a1c_mean_quant$',i,'[ups_13_feat], ups_13_feat),
    as.data.frame(v1c_mean_quant$',i,'[ups_13_feat], ups_13_feat),
    as.data.frame(m1_mean_quant$',i,'[ups_13_feat], ups_13_feat),
    as.data.frame(s1_mean_quant$',i,'[ups_13_feat], ups_13_feat),
    as.data.frame(cgg_mean_quant$',i,'[ups_13_feat], ups_13_feat)
  )')
  
  eval(parse(text=command))
    
  hist_in <- apply(hist_dat, 1, function(x) mean(x, na.rm = TRUE))
  hist(hist_in, breaks = 13, main = i)
}


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