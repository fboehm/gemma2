#' Calculate eigendecomposition and return ordered eigenvalues and eigenvectors
#'
#' @param spd a semi-positive definite matrix
#' @param decreasing argument passed to order()
#' @export
eigen2 <- function(spd, decreasing = FALSE){
  eigen(spd) -> foo
  bar <- foo
  bar$values <- foo$values[order(foo$values, decreasing = decreasing)]
  bar$vectors <- foo$vectors[, order(foo$values, decreasing = decreasing)]
  return(bar)
}
