library(gemma2)
library(readr)
#read_tsv("~/GEMMA/example/mouse100.pheno.txt", col_names = FALSE) -> pheno
#read_tsv("~/GEMMA/example/output/mouse100.cXX.txt", col_names = FALSE)[, 1:100] -> kinship
readr::read_tsv(system.file("extdata", "mouse100.pheno.txt", package = "gemma2"), col_names = FALSE) -> pheno
readr::read_tsv(system.file("extdata", "mouse100.cXX.txt", package = "gemma2"), col_names = FALSE)[, 1:100] -> kinship

e_out <- eigen2(kinship)
as.matrix(pheno[, c(1, 6)]) -> phe16

MphEM(max_iter = 10, eval = e_out$values, X = matrix(c(-10, rep(0, 99)), nrow = 1),
      Y = t(phe16) %*% e_out$vectors,
      V_g = matrix(c(1.91352, 0, 0, 0.530827), nrow = 2),
      V_e = matrix(c( 0.320028, 0, 0, 0.561589), nrow = 2)
)-> foo

context("Testing MphEM() with two phenotypes and intercept only (no genetic data)")



test_that("First calculations of for Vg and for Ve are accurate", {
  expect_equal(foo[[1]][[2]], matrix(c(1.91352, 0.0700492, 0.0700492, 0.530827), nrow = 2, ncol = 2, byrow = FALSE), tolerance = 0.000001)
  expect_equal(foo[[1]][[3]], matrix(c(0.320028, 0.0691272, 0.0691272, 0.561589), nrow = 2, ncol = 2, byrow = FALSE), tolerance = 0.000001)
})
