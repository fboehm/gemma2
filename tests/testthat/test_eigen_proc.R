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
