#' Update U matrix
#'
#' @param OmegaE the OmegaE matrix, calculated in calc_omega
#' @param UltVehiY matrix
#' @param UltVehiBX matrix
#' @family expectation-maximization functions
#' @examples
#' readr::read_tsv(system.file("extdata",
#' "mouse100.pheno.txt",
#' package = "gemma2"),
#' col_names = FALSE) -> pheno
#' phe16 <- as.matrix(pheno[, c(1, 6)])
#' as.matrix(readr::read_tsv(system.file("extdata",
#' "mouse100.cXX.txt",
#' package = "gemma2"),
#' col_names = FALSE)[, 1:100]) -> kinship
#' eigen2(kinship) -> e2_out
#' e2_out$values -> eval
#' e2_out$vectors -> U
#' eigen_proc(V_g = diag(c(1.91352, 0.530827)),
#' V_e = diag(c(0.320028, 0.561589))) -> ep_out
#' UltVehi <- ep_out[[3]]
#' calc_omega(eval, ep_out$D_l) -> co_out
#' update_u(OmegaE = co_out[[2]],
#'         UltVehiY = UltVehi %*% t(phe16),
#'         UltVehiBX = matrix(c(-0.71342, -0.824482),
#'         ncol = 1) %*% t(rep(1, 100))
#' )
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
#' @family expectation-maximization functions
update_e <- function(UltVehiY, UltVehiBX, UltVehiU){
  UltVehiE <- UltVehiY - UltVehiBX - UltVehiU
  return(UltVehiE)
}

#' Update B for restricted log likelihood
#'
#' @param xHiy vector
#' @param Qi Q inverse matrix
#' @param d_size number of traits
#' @family expectation-maximization functions
UpdateRL_B <- function(xHiy, Qi, d_size){
  nrow(Qi) -> dc_size
  c_size <- dc_size / d_size
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
#' @param tol a positive number indicating tolerance to be passed to isSymmetric()
#' @family expectation-maximization functions
update_v <- function(eval, U, E, Sigma_uu, Sigma_ee, tol = 1 / 10000){
  stopifnot(isSymmetric(Sigma_uu, tol = tol), isSymmetric(Sigma_ee, tol = tol))
  n_size <- length(eval)
  d_size <- nrow(U)
  V_g <- matrix(0, nrow = d_size, ncol = d_size)
  V_e <- V_g
  for (k in 1:n_size){
    delta <- eval[k]
    #if (delta != 0){
      U_col <- U[, k]
      V_g <- V_g + U_col %*% t(U_col) / delta
    #}
  }
  V_e <- E %*% t(E)
  V_g <- V_g + Sigma_uu
  V_e <- V_e + Sigma_ee
  V_g <- V_g / n_size
  V_e <- V_e / n_size
  return(list(V_g, V_e))
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
calc_sigma <- function(eval, D_l, X, OmegaU, OmegaE, UltVeh, Qi){
  n_size <- length(eval)
  c_size <- nrow(X)
  d_size <- length(D_l)
  dc_size <- nrow(Qi)
  Sigma_ee <- matrix(0, nrow = d_size, ncol = d_size)
  Sigma_uu <- Sigma_ee
  for (k in 1:n_size){
    OmegaU_col <- OmegaU[, k]
    OmegaE_col <- OmegaE[, k]
    diag(Sigma_uu) <- diag(Sigma_uu) + OmegaU_col
    diag(Sigma_ee) <- diag(Sigma_ee) + OmegaE_col
  }
  M_u <- matrix(0, nrow = dc_size, ncol = d_size)
  M_e <- M_u
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
  } # end loop in k
  M <- Sigma_uu %*% UltVeh
  Sigma_uu <- t(UltVeh) %*% M
  M <- Sigma_ee %*% UltVeh
  Sigma_ee <- t(UltVeh) %*% M
  return(list(Sigma_ee, Sigma_uu))
}

