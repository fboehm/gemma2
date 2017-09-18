#' Calculate XHiY for two loci
#'
#' @param eval vector of eigenvalues from the decomposition of the relatedness matrix
#' @param D_l vector of length d_size
#' @param X1 design matrix for locus 1
#' @param X2 design matrix for locus 2
#' @param UltVehiY a matrix
#' @export

calc_XHiY_2loci <- function (eval, D_l, X1, X2, UltVehiY)
{
  stopifnot(length(eval) == ncol(X1), is.vector(eval), is.vector(D_l))
  n_size <- length(eval)
  c_size <- nrow(X1)
  d_size <- length(D_l)
  X <- array(data = c(X1, X2),
             dim = c(nrow(X1), ncol(X1), 2))
  xHiy <- rep(0, d_size * c_size)
  for (i in 1:d_size) {
    dl <- D_l[i]
    for (j in 1:c_size) {
      d <- 0
      for (k in 1:n_size) {
        x <- X[j, k, i] # NOTE THREE INDICES HERE
        y <- UltVehiY[i, k]
        delta <- eval[k]
        d <- d + x * y/(delta * dl + 1)
      }
      xHiy[(j - 1) * d_size + i] <- d
    }
  }
  return(xHiy)
}
