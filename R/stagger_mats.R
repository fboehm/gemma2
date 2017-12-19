#' Stagger matrices of same size within a larger, block-diagonal matrix
#'
#' @param X1 first matrix
#' @param X2 second matrix, of same size as first
#' @export

stagger_mats <- function(...){
  out <- as.matrix(Matrix::bdiag(...))
  return(out)
}
