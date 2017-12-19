#' Calculate Sigma_ee and Sigma_uu matrices for 2 loci
#'
#' @param eval eigenvalues vector from decomposition of relatedness matrix
#' @param D_l vector
#' @param X1 design matrix
#' @param X2 design matrix
#' @param OmegaU matrix
#' @param OmegaE matrix
#' @param UltVeh matrix
#' @param Qi inverse of Q matrix
#' @export
calc_sigma_2loci <- function (eval, D_l, X1, X2, OmegaU, OmegaE, UltVeh, Qi)
{
  n_size <- length(eval)
  c_size <- nrow(X1)
  d_size <- length(D_l)
  dc_size <- nrow(Qi)
  Sigma_ee <- matrix(0, nrow = d_size, ncol = d_size)
  Sigma_uu <- Sigma_ee
  for (k in 1:n_size) {
    OmegaU_col <- OmegaU[, k]
    OmegaE_col <- OmegaE[, k]
    diag(Sigma_uu) <- diag(Sigma_uu) + OmegaU_col
    diag(Sigma_ee) <- diag(Sigma_ee) + OmegaE_col
  }
  M_u <- matrix(0, nrow = dc_size, ncol = d_size)
  M_e <- M_u
  X <- array(c(X1, X2), dim = c(nrow(X1), ncol(X1), 2))
  for (k in 1:n_size) {
    delta <- eval[k]
    for (i in 1:d_size) {
      dl <- D_l[i]
      for (j in 1:c_size) {
        x <- X[j, k, i]
        d <- x/(delta * dl + 1)
        M_e[(j - 1) * d_size + i, i] <- d
        M_u[(j - 1) * d_size + i, i] <- d * dl
      }
    }
    QiM <- Qi %*% M_u
    Sigma_uu <- Sigma_uu + t(M_u) %*% QiM * delta
    QiM <- Qi %*% M_e
    Sigma_ee <- Sigma_ee + t(M_e) %*% QiM
  }
  M <- Sigma_uu %*% UltVeh
  Sigma_uu <- t(UltVeh) %*% M
  M <- Sigma_ee %*% UltVeh
  Sigma_ee <- t(UltVeh) %*% M
  return(list(Sigma_ee, Sigma_uu))
}

