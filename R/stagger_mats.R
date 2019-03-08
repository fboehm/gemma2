#' Stagger matrices within a larger, block-diagonal matrix
#'
#' @param ... one or more matrices, separated by commas
#' @return a block-diagonal matrix, with the inputted matrices as blocks on the diagonal.
#' @examples
#' foo <- matrix(rnorm(40000), ncol = 8)
#' block_diag <- stagger_mats(foo, foo)
#' dim(foo)
#' dim(block_diag)
#' @export



stagger_mats <- function(...){
  return(as.matrix(Matrix::bdiag(...)))
}
