library(gemma2)
library(readr)
read_tsv("~/GEMMA/example/mouse100.pheno.txt", col_names = FALSE) -> pheno
read_tsv("~/GEMMA/example/output/mouse100.cXX.txt", col_names = FALSE)[, 1:100] -> kinship
e_out <- eigen2(kinship)
library(tidyverse)
as.matrix(pheno[, c(1, 6)]) -> phe16

MphEM(max_iter = 10, eval = e_out$values, X = matrix(c(-10, rep(0, 99)), nrow = 1),
      Y = t(phe16) %*% e_out$vectors,
      V_g = matrix(c(1.91352, 0, 0, 0.530827), nrow = 2),
      V_e = matrix(c( 0.320028, 0, 0, 0.561589), nrow = 2)
)-> foo

context("Testing MphEM() with two phenotypes and intercept only (no genetic data)")



test_that("First calculation of off-diagonal term for Vg is accurate", {
  expect_equal(foo[[1]][[2]][1, 2], 0.07004918)

})
