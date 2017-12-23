#' Center a relatedness matrix, after Zhou's GEMMA function CenterMatrix
#'
#' @param mat a relatedness matrix
#' @export

center_kinship <- function(mat){
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
#// Center the matrix G.
#void CenterMatrix(gsl_matrix *G) {
#  double d;
#  gsl_vector *w = gsl_vector_alloc(G->size1);
#  gsl_vector *Gw = gsl_vector_alloc(G->size1);

#
#
#
#
#  for (size_t i = 0; i < G->size1; ++i) {
#    for (size_t j = 0; j < i; ++j) {
#      d = gsl_matrix_get(G, j, i);
#      gsl_matrix_set(G, i, j, d);
#    }
#  }#
#
#
#  gsl_vector_free(w);
#  gsl_vector_free(Gw);#
#
#
#  return;
#}