#' Calculate Qi (inverse of Q) and log determinant of Q for 2 loci (which might be the same locus)
#'
#' @param eval vector of eigenvalues from decomposition of relatedness matrix
#' @param D_l eigenvalues vector of length d_size
#' @param X a 16 by n design matrix with genotype probabilities for two loci
#' @export
# //Qi=(\sum_{k=1}^n x_kx_k^T\otimes(delta_k*Dl+I)^{-1} )^{-1} according to Zhou's mvlmm.cpp file
calc_qi2 <- function(eval, D_l, X){
  n_size <- ncol(X)
  I2 <- diag(1, 2)
  Q <- 0
  for (k in 1:n_size){
    xk <- X[, k]
    delta_k <- eval[k]
    Q <- Q + xk %*% t(xk) %x% solve(delta_k * D_l + I2)
  }
  Qi <- solve(Q)
  detQ <- det(Q)
  lndetQ <- log(detQ)
  return(list(Qi, lndetQ))
}
