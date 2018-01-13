library(gemma2)
context("Testing calc_sigma")

readr::read_tsv(system.file("extdata", "mouse100.pheno.txt", package = "gemma2"), col_names = FALSE) -> pheno
phe16 <- as.matrix(pheno[, c(1, 6)])
as.matrix(readr::read_tsv(system.file("extdata", "mouse100.cXX.txt", package = "gemma2"), col_names = FALSE)[, 1:100]) -> kinship

eigen2(kinship) -> eout
eout$values -> eval
eout$vectors -> U

V_g <- diag(c(1.91352, 0.530827))
V_e <- diag(c(0.320028, 0.561589))

X <- t(rep(1, 100)) %*% U

ep_out <- eigen_proc(V_g, V_e)
ep_out[[1]] -> logdet_Ve
ep_out[[2]] -> UltVeh
ep_out[[3]] -> UltVehi
ep_out[[4]] -> D_l

cq_out <- calc_qi(eval, D_l, X)
cq_out[[1]] -> Qi

co_out <- calc_omega(eval, D_l)
co_out[[1]] -> OmegaU
co_out[[2]] -> OmegaE

calc_sigma(eval, D_l, X, OmegaU, OmegaE, UltVeh, Qi) -> cs_out
cs_out[[1]] -> Sigma_ee
cs_out[[2]] -> Sigma_uu

test_that("Results of gemma2 equal those of GEMMA v 0.97", {
  expect_equal(Sigma_ee, diag(c(18.559, 12.3672)), tolerance = 0.0001)
  expect_equal(Sigma_uu, diag(c(82.2973, 41.9238)), tolerance = 0.0001)
})

