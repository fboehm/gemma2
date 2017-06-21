#' Update U matrix
#'
#' @param OmegaE the OmegaE matrix, calculated in calc_omega
#' @param UltVehiY matrix
#' @param UltVehiBX matrix
#' @export
update_u <- function(OmegaE, UltVehiY, UltVehiBX){
  UltVehiU <- UltVehiY
  UltVehiU <- UltVehiU - UltVehiBX
  UltVehiU <- UltVehiU * OmegaE
  return(UltVehiU)
}

#' Update E
#'
#' @param UltVehiY matrix of transformed Y values
#' @param UltVehiBX matrix of transformed BX values
#' @param UltVehiU matrix of transformed U values
#' @export
update_e <- function(UltVehiY, UltVehiBX, UltVehiU){
  UltVehiY -> UltVehiE
  UltVehiE <- UltVehiE - UltVehiBX
  UltVehiE <- UltVehiE - UltVehiU
  return(UltVehiE)
}

#' Update restricted log likelihood
#'
#' @param xHiy vector
#' @param Qi Q inverse matrix
#' @param d_size number of traits
UpdateRL_B <- function(xHiy, Qi, d_size = 2){
  nrow(Qi) -> dc_size
  c_size <- dc_size / d_size
  b <- vector(length = dc_size)
  b <- Qi %*% xHiy
  UltVehiB <- matrix(nrow = d_size, ncol = c_size)
  for (i in 1:c_size){
    b_subcol <- b[(1 + (i - 1) * d_size):(i * d_size)]
    b_subcol -> UltVehiB[, i] # could use as.matrix here
  }
  return(UltVehiB)
}

#' Update V_e and V_g
#'
#' @param eval vector of eigenvalues from eigendecomposition of relatedness matrix
#' @param U matrix
#' @param E matrix
#' @param Sigma_uu matrix
#' @param Sigma_ee matrix
#' @export
update_v <- function(eval, U, E, Sigma_uu, Sigma_ee){
  n_size <- length(eval)
  d_size <- nrow(U)
  V_g <- matrix(0, nrow = d_size, ncol = d_size)
  V_e <- V_g
  for (k in 1:n_size){
    delta <- eval[k]
    if (delta != 0){
      U_col <- U[, k]
      V_g <- V_g + U_col %*% t(U_col) / delta
    }
  }
  V_e <- E %*% t(E)
  V_g <- V_g + Sigma_uu
  V_e <- V_e + Sigma_ee
  V_g <- V_g / n_size
  V_e <- V_e / n_size
  return(list(V_e, V_g))
}

#' Calculate Sigma_ee and Sigma_uu matrices
#'
#' @param eval eigenvalues vector from decomposition of relatedness matrix
#' @param D_l vector
#' @param X design matrix
#' @param OmegaU matrix
#' @param OmegaE matrix
#' @param UltVeh matrix
#' @param Qi inverse of Q matrix
#' @param func_name indicator for restricted or unrestricted likelihood
#' @export
calc_sigma <- function(eval, D_l, X, OmegaU, OmegaE, UltVeh, Qi, func_name = "R"){
  n_size <- length(eval)
  c_size <- nrow(X)
  d_size <- length(D_l)
  dc_size <- nrow(Qi)
  Sigma_ee <- matrix(0, nrow = d_size, ncol = d_size)
  Sigma_uu <- Sigma_ee
  for (k in 1:n_size){
    OmegaU_col <- OmegaU[, k]
    OmegaE_col <- OmegaE[, k]
    diag(Sigma_ee) <- diag(Sigma_ee) + OmegaE_col
    diag(Sigma_uu) <- diag(Sigma_uu) + OmegaU_col
  }
  if (func_name == "R"){
    M_u <- matrix(0, nrow = dc_size, ncol = d_size)
    M_e <- M_u
    QiM <- M_e
    for (k in 1:n_size){
      delta <- eval[k]
      for (i in 1:d_size){
        dl <- D_l[i]
        for (j in 1:c_size){
          x <- X[j, k]
          d <- x / (delta * dl + 1)
          M_e[(j - 1) * d_size + i, i] <- d
          M_u[(j - 1) * d_size + i, i] <- d * dl

        }
      }
      QiM <- Qi %*% M_u
      Sigma_uu <- Sigma_uu + t(M_u) %*% QiM * delta
      QiM <- Qi %*% M_e
      Sigma_ee <- Sigma_ee + t(M_e) %*% QiM
    }
  }
  M <- Sigma_uu %*% UltVeh
  Sigma_uu <- t(UltVeh) %*% M
  M <- Sigma_ee %*% UltVeh
  Sigma_ee <- t(UltVeh) %*% M
  return(list(Sigma_ee, Sigma_uu))
}

