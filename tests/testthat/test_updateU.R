library(gemma2)
context("testing update_u")

readr::read_tsv(system.file("extdata", "mouse100.pheno.txt", package = "gemma2"), col_names = FALSE) -> pheno
phe16 <- as.matrix(pheno[, c(1, 6)])

as.matrix(readr::read_tsv(system.file("extdata", "mouse100.cXX.txt", package = "gemma2"), col_names = FALSE)[, 1:100]) -> kinship
eigen2(kinship) -> e2_out
e2_out$values -> eval
e2_out$vectors -> U
eigen_proc(V_g = diag(c(1.91352, 0.530827)), V_e = diag(c(0.320028, 0.561589))) -> ep_out

UltVehi <- ep_out[[3]]

calc_omega(eval, ep_out$D_l) -> co_out
update_u(OmegaE = co_out[[2]],
         UltVehiY = UltVehi %*% t(phe16),
         UltVehiBX = matrix(c(-0.71342, -0.824482), ncol = 1) %*% t(rep(1, 100))
         ) -> UltVehiU

test_that("update_u output matches GEMMAv0.97 output for intercept-only model", {
  expect_equal(UltVehiU[, 1], c(0, 0), tolerance = 0.001)
  expect_length(UltVehiU, 200)
})
