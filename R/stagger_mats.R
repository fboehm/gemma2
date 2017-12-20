#' Stagger matrices within a larger, block-diagonal matrix
#'
#' @param ... one or more matrices, separated by commas
#' @export

stagger_mats <- function(...){
  out <- as.matrix(Matrix::bdiag(...))
  return(out)
}
