#' Center a kinship matrix before eigendecomposition
#'
#' @param mat a kinship matrix
#' @export

center_matrix <- function(mat){
  #  gsl_vector_set_all(w, 1.0);
  w <- rep(1, nrow(mat))
  #  gsl_blas_dgemv(CblasNoTrans, 1.0, G, w, 0.0, Gw);
  mat %*% w -> matw
  #  gsl_blas_dsyr2(CblasUpper, -1.0 / (double)G->size1, Gw, w, G);
  foo <- mat - (matw %*% t(w) + w %*% t(matw)) / nrow(mat)
  #  gsl_blas_ddot(w, Gw, &d);
  d <- t(w) %*% matw
  #  gsl_blas_dsyr(CblasUpper, d / ((double)G->size1 * (double)G->size1), w, #G);
  out <- foo + (w %*% t(w)) * as.numeric(d / (nrow(mat)^2))
  return(out)
}
