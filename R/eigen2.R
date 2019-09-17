#' Calculate eigendecomposition and return ordered eigenvalues and eigenvectors
#'
#' @param spd a semi-positive definite matrix
#' @param decreasing argument passed to order()
#' @export
#' @examples
#' readr::read_tsv(system.file("extdata", "mouse100.cXX.txt", package = "gemma2"), col_names = FALSE)[, 1:100] -> kinship
#' e_out <- eigen2(as.matrix(kinship))
#' @return a list with 2 components, the eigenvalues and the eigenvectors
eigen2 <- function(spd, decreasing = FALSE){
  eigen(spd) -> foo
  bar <- foo
  bar$values <- foo$values[order(foo$values, decreasing = decreasing)]
  bar$vectors <- foo$vectors[, order(foo$values, decreasing = decreasing)]
  return(bar)
}
