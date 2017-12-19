#' Calculate log likelihood
#'
#' @param eval eigenvalues vector from decomposition of relatedness matrix
#' @param D_l vector of eigenvalues from decomposition of Ve matrix
#' @param Qi inverse of Q matrix
#' @param UltVehiY matrix of (transformed) Y values
#' @param xHiy vector
#' @export
MphCalcLogL <- function(eval, D_l, Qi, UltVehiY, xHiy){
  n_size <- length(eval)
  d_size <- length(D_l) # d is number of phenotypes
  dc_size <- nrow(Qi)
  logl <- 0
  for (k in 1:n_size){
    delta <- eval[k]
    for (i in 1:d_size){
      y <- UltVehiY[i, k]
      dl <- D_l[i]
      d <- delta * dl + 1
      logl <- logl + y^2 / d + log(d)
    }
  }
  Qiv <- Qi %*% xHiy
  d <- t(xHiy) %*% Qiv
  stopifnot(length(d) == 1)
  logl <- logl - d
  return(- 0.5 * logl)
}

