#' Calculate Qi (inverse of Q) and log determinant of Q for two loci
#'
#' @param eval vector of eigenvalues from decomposition of relatedness matrix
#' @param D_l vector of length d_size
#' @param X1 a design matrix
#' @param X2 a design matrix
#' @export
calc_qi_2loci <- function (eval, D_l, X1, X2)
{
  n_size <- length(eval)
  d_size <- length(D_l)
  c_size <- nrow(X1)
  dc_size <- d_size * c_size
  # define X
  X <- array(c(X1, X2), dim = c(nrow(X1), ncol(X1), 2))
  Q <- matrix(0, nrow = dc_size, ncol = dc_size)
  for (i in 1:c_size) {
    for (j in 1:c_size) {
      for (l in 1:d_size) {
        dl <- D_l[l]
        #if (j < i) {
        #  d <- Q[(j - 1) * d_size + l, (i - 1) * d_size +
        #    l]
        #}
        #else {
        d <- 0
        for (k in 1:n_size) {
          d1 <- X[i, k, l]
          d2 <- X[j, k, l]
          delta <- eval[k]
          d <- d + d1 * d2/(dl * delta + 1)
        }
        #}
        Q[(i - 1) * d_size + l, (j - 1) * d_size + l] <- d
      }
    }
  }
  Qi <- solve(Q)
  detQ <- det(Q)
  lndetQ <- log(detQ)
  return(list(Qi, lndetQ))
}
