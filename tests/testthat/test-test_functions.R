source("/home/nperumal/AllenBrainSC/scripts/functions.R")
library(testthat)
library(parallel)


meta_test <- data.frame('sample_name' = c('sample1','sample2','sample3','sample4','sample5'),
                        'region_label' = c('MTG','MTG','MTG','V1C','V1C'),
                        'broad_cell' = c('Exc','Inh','Unk','Exc','Inh'),
                        'cell_type_alias_label' = c('Exc_L1','Inh_L1','Unk_L1','Exc_L2','Inh_L2')
)

gene_exp_test <- data.frame('sample1' = 1:20,
                            'sample2' = 21:40,
                            'sample3' = 41:60,
                            'sample4' = 61:80,
                            'sample5' = 71:90)
rownames(gene_exp_test) <- c('gene1','gene2','gene3','gene4','gene5','gene6','gene7',
                             'gene8','gene9','gene10','gene11','gene12','gene13','gene14',
                             'gene15','gene16','gene17','gene18','gene19','gene20')




# CPM calculation tests
test_that("CPM calculation works", {
  expect_equal(CPM(c(15,0)), c(1000000,0)) # Checks for value
  suppressWarnings(expect_setequal(CPM(as.data.frame(c(15,0))),
                                   as.data.frame(c(1000000,0)))) # Checks output type
})


# Testing the new quant pruning function
# Check for numbers of items to see if the function work
# Need to fix proper pruning
test_that("Error messages missing parameters",{
  expect_error(prune_quant(met = meta, 
                          region = 'MTG',
                          is.specific = FALSE, 
                          result = 'feature_name',
                          cell_pop_size = 50),
               'Must specify expression matrix')
  
  expect_error(prune_quant(exp = gene_exp, 
                          region = 'MTG',
                          is.specific = FALSE, 
                          result = 'feature_name',
                          cell_pop_size = 50),
               'Must specify meta dataframe')
  
  expect_error(prune_quant(exp = gene_exp,
                              met = meta, 
                              is.specific = FALSE, 
                              result = 'feature_name',
                           cell_pop_size = 50),
               'Must specify a region')
  
  expect_error(prune_quant(exp = gene_exp,
                          met = meta,
                          region = 'MTG',
                          result = 'feature_name',
                          cell_pop_size = 50),
               'Must specify is.specific with boolean value')
  
  expect_error(prune_quant(exp = gene_exp,
                          met = meta,
                          region = 'MTG',
                          is.specific = FALSE, 
                          result = 'Random_input',
                          cell_pop_size = 50),
               'Result must be either feature_name or exp_matrix')
  
})


test_that('Confirming prune_quant output with test data',{
  
  expect_true(length(prune_quant(exp = gene_exp_test,
                           met = meta_test,
                           region = 'MTG',
                           is.specific = TRUE,
                           result = 'feature_name',
                           cell_pop_size = 1)) == 3)
  
  
  expect_true(length(prune_quant(exp = gene_exp_test,
                                 met = meta_test,
                                 region = 'MTG',
                                 is.specific = TRUE,
                                 result = 'feature_name',
                                 cell_pop_size = 1)$Exc) == 1)

})



