library(gemma2)
context("Testing calc_XHiY")

readr::read_tsv(system.file("extdata", "mouse100.pheno.txt", package = "gemma2"), col_names = FALSE) -> pheno
phe16 <- as.matrix(pheno[, c(1, 6)])
as.matrix(readr::read_tsv(system.file("extdata", "mouse100.cXX.txt", package = "gemma2"), col_names = FALSE)[, 1:100]) -> kinship
library(gemma2)
eigen2(kinship) -> eout
eout$values -> eval
eout$vectors -> U


UltVehi <- matrix(c(0, -1.76769, -1.334414, 0), nrow = 2, byrow = FALSE) # from output of eigen_proc()

calc_XHiY(eval = eval, D_l = c(0.9452233, 5.9792268),
          X = rep(1, 100) %*% U,
          UltVehiY = UltVehi %*% t(phe16) %*% U) -> foo

test_that("calc_XHiY outputs match those of GEMMAv0.97", {
  expect_equal(foo, c(-71.342, -82.4482), tolerance = 0.0001)
})
