library(gemma2)
context("Testing eigen_proc")
test_that("outputted logdetVe value is length 1", {
  expect_equal(length(eigen_proc(matrix(c(1, 0, 0, 1), nrow = 2),
                                 matrix(c(1, 0, 0, 1), nrow = 2))[[1]]),
               1)
})

test_that("outputted logdetVe value is non-negative", {
  expect_gte(eigen_proc(matrix(c(1, 0, 0, 1), nrow = 2),
                                 matrix(c(1, 0, 0, 1), nrow = 2))[[1]], 0)
  expect_gt(eigen_proc(matrix(c(5, 0, 0, 1), nrow = 2),
                        matrix(c(5, 0, 0, 1), nrow = 2))[[1]], 0)
})

test_that("second and third outputs have same dimensions", {
  expect_equal(nrow(eigen_proc(matrix(c(1, 0, 0, 1), nrow = 2),
             matrix(c(1, 0, 0, 1), nrow = 2))[[2]]),
             nrow(eigen_proc(matrix(c(1, 0, 0, 1), nrow = 2),
                             matrix(c(1, 0, 0, 1), nrow = 2))[[3]]))
  expect_equal(ncol(eigen_proc(matrix(c(1, 0, 0, 1), nrow = 2),
                               matrix(c(1, 0, 0, 1), nrow = 2))[[2]]),
               ncol(eigen_proc(matrix(c(1, 0, 0, 1), nrow = 2),
                               matrix(c(1, 0, 0, 1), nrow = 2))[[3]]))
})

test_that("fourth output has proper length", {
  expect_length(eigen_proc(matrix(c(1, 0, 0, 1), nrow = 2),
                           matrix(c(1, 0, 0, 1), nrow = 2))[[4]], 2)
})

test_that("entries of fourth output are all non-negative", {
  expect_gte(eigen_proc(matrix(c(1, 0, 0, 1), nrow = 2),
                           matrix(c(1, 0, 0, 1), nrow = 2))[[4]][1], 0)
  expect_gte(eigen_proc(matrix(c(1, 0, 0, 1), nrow = 2),
                        matrix(c(1, 0, 0, 1), nrow = 2))[[4]][2], 0)
})
foo <- abs(rnorm(4))
test_that("first entry of output is correctly calculated", {
  expect_equal(eigen_proc(matrix(c(1, 0, 0, 1), nrow = 2),
                          matrix(c(1, 0, 0, 1), nrow = 2))[[1]], 0)
  expect_equal(eigen_proc( matrix(c(foo[1], 0, 0, foo[2]), nrow = 2),
                          matrix(c(foo[3],0, 0, foo[4]), nrow = 2))[[1]],
               log(det(matrix(c(foo[3],0, 0, foo[4]), nrow = 2))))
})

# Compare with GEMMAv0.97 outputs from EigenProc
test_that("simple case agrees with output from GEMMA's EigenProc",{
  expect_equal(eigen_proc(V_g = diag(c(1.91352, 0.530827)), V_e = diag(c(0.320028, 0.561589)))[[1]], -1.716332,
               tolerance = 0.000001)
  expect_equal(eigen_proc(V_g = diag(c(1.91352, 0.530827)), V_e = diag(c(0.320028, 0.561589)))[[2]],
               matrix(c(0, -0.5657102, -0.7493924, 0), nrow = 2, byrow = FALSE),
               tolerance = 0.000001
               )
  expect_equal(eigen_proc(V_g = diag(c(1.91352, 0.530827)), V_e = diag(c(0.320028, 0.561589)))[[3]],
               matrix(c(0, -1.76769, -1.334414, 0), nrow = 2, byrow = FALSE),
               tolerance = 0.000001
               )
  expect_equal(eigen_proc(V_g = diag(c(1.91352, 0.530827)), V_e = diag(c(0.320028, 0.561589)))[[4]],
               c(0.9452233, 5.9792268),
               tolerance = 0.000001
               )
})
