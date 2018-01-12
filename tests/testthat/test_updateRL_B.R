library(gemma2)
context("Testing updateRL_B")

UpdateRL_B(xHiy = c(-71.342, -82.4482), Qi = diag(rep(0.01, 2)), d_size = 2) -> bar

test_that("updateRL_B output matches that of GEMMA v 0.97", {
  expect_equal(as.vector(bar), c(-0.71342, -0.824482), tolerance = 0.0001)
})
