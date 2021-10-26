source("/home/nperumal/AllenBrainSC/scripts/functions.R")
library(testthat)

# CPM calculation tests
test_that("CPM calculation works", {
  expect_equal(CPM(c(15,0)), c(1000000,0)) # Checks for value
  suppressWarnings(expect_setequal(CPM(as.data.frame(c(15,0))),
                                   as.data.frame(c(1000000,0)))) # Checks output type
})


# Testing the new quant pruning function
test_that("Error messages missing parameters",{
  expect_error(test_quant(met = meta, 
                          region = 'MTG',
                          is.specific = FALSE, 
                          result = 'feature_name'),
               'Must specify expression matrix')
  
  expect_error(test_quant(exp = gene_exp, 
                          region = 'MTG',
                          is.specific = FALSE, 
                          result = 'feature_name'),
               'Must specify meta dataframe')
  
  expect_error(test_quant(exp = gene_exp,
                              met = meta, 
                              is.specific = FALSE, 
                              result = 'feature_name'),
               'Must specify a region')
  
  expect_error(test_quant(exp = gene_exp,
                          met = meta,
                          region = 'MTG',
                          result = 'feature_name'),
               'Must specify is.specific with boolean value')
  
  expect_error(test_quant(exp = gene_exp,
                          met = meta,
                          region = 'MTG',
                          is.specific = FALSE, 
                          result = 'Random_input'),
               'Result must be either feature_name or exp_matrix')
  
})

test_that('Output values are the same',{
  expect_identical(test_quant(exp = gene_exp,
                             met = meta,
                             region = 'MTG',
                             is.specific = FALSE,
                             result = 'feature_name')$Exc,
                  
                  rownames(filter_dat(exp = gene_exp, met = meta, 
                                      pcnt = 0.75, value = 1, 
                                      region = 'MTG', is.broad = TRUE, 
                                      cell_t = 'Exc')$exp)
               )
})

