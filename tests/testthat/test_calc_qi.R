library(gemma2)
context("Testing calc_qi")

as.matrix(readr::read_tsv(system.file("extdata", "mouse100.cXX.txt", package = "gemma2"), col_names = FALSE)[, 1:100]) -> kinship
eigen2(kinship) -> e2_out
e2_out$values -> eval
e2_out$vectors -> U
eigen_proc(V_g = diag(c(1.91352, 0.530827)), V_e = diag(c(0.320028, 0.561589))) -> ep_out

calc_qi(eval = eval, D_l = ep_out[[4]], X = t(rep(1, 100)) %*% U) -> cq_out

test_that("logdetVe and Qi match that from GEMMAv0.97 for intercept-only model",{
          expect_equal(cq_out[[1]], diag(rep(0.01, 2)), tolerance = 0.0001)
          expect_equal(cq_out[[2]], 9.21034, tolerance = 0.00001)
          })
