#' Calculate XHiY
#'
#' @param eval vector of eigenvalues from the decomposition of the relatedness matrix
#' @param D_l vector of length d_size
#' @param X design matrix
#' @param UltVehiY a matrix
#' @export
calc_XHiY <- function(eval, D_l, X, UltVehiY){
  n_size <- length(eval)
  c_size <- nrow(X)
  d_size <- length(D_l)
  xHiy <- numeric(length = d_size * c_size)
  for (i in 1:d_size){
    dl <- D_l[i]
    for (j in 1:c_size){
      d <- 0
      for (k in 1:n_size){
        x <- X[j, k]
        y <- UltVehiY[i, k]
        delta <- eval[k]
        d<- d + x * y / (delta * dl + 1)
      }
    }
    xHiy[(j-1) * d_size + i] <- d
  }
  return(xHiy)
}

#' Eigendecomposition procedure for Vg and Ve
#'
#' @param Vg a d_size by d_size covariance matrix
#' @param Ve a d_size by d_size covariance matrix
#' @export
eigen_proc <- function(Vg, Ve){
  eigenVe <- eigen(Ve)
  eigenVe$values -> Dl
  eigenVe$vectors -> Ul
  logDl <- log(Dl)
  sum(logDl) -> logdet_Ve
  V_e_h <- matrix(0, nrow = length(Dl), ncol = length(Dl))
  V_e_hi <- V_e_h
  for (i in 1:length(Dl)){
    V_e_h <- V_e_h + sqrt(Dl)[i] * Ul[, i] %*% t(Ul[, i])
    V_e_hi <- V_e_hi + Ul[, i] %*% t(Ul[, i]) /  sqrt(Dl)[i]
  }
  VgVehi <- Vg %*% V_e_hi
  Lambda <- V_e_hi %*% VgVehi
  eigen(Lambda) -> eigenLambda
  Ul <- eigenLambda$vectors
  Dl <- eigenLambda$values
  Dl[Dl < 0] <- 0
  UltVeh <- t(Ul) %*% V_e_h
  UltVehi <- t(Ul) %*% V_e_hi
  return(list(logdet_Ve, UltVeh, UltVehi, Dl))
}


#' Calculate Qi (inverse of Q) and log determinant of Q
#'
#' @param eval vector of eigenvalues from decomposition of relatedness matrix
#' @param D_l vector of length d_size
#' @param X a design matrix
#' @export
calc_qi <- function(eval, D_l, X){
  n_size <- length(eval)
  d_size <- length(D_l)
  c_size <- nrow(X)# what is c_size? it's the number of rows in the transposed genotypes matrix
  # transposed genotypes matrix is c by n
  dc_size <- d_size * c_size
  Q <- matrix(0, nrow = dc_size, ncol = dc_size)

  for (i in 1:c_size){
    for (j in 1:c_size){
      for (l in 1:d_size){
        dl <- D_l[l]
        if (j < i){
          d <- Q[(j - 1) * d_size + l, (i - 1) * d_size + l]
        } else {
          d <- 0
          for (k in 1:n_size){
            d1 <- X[i, k]
            d2 <- X[j, k]
            delta <- eval[k]
            d <- d + d1 * d2 / (dl * delta + 1)
          }
        }
        Q[(i -1) * d_size + l, (j - 1) * d_size + l] <- d
      }
    }
  }
  Qi <- solve(Q)
  detQ <- det(Q)
  lndetQ <- log(detQ)
  return(list(Qi, lndetQ))
}

#' Calculate Omega matrices
#'
#' @param eval vector of eigenvalues from decomposition of relatedness matrix
#' @param D_l vector of length d_size
#' @export
calc_omega <- function(eval, D_l){
  n_size <- length(eval)
  d_size <- length(D_l)
  OmegaU <- matrix(nrow = d_size, ncol = n_size)
  OmegaE <- OmegaU
  for (k in 1:n_size){
    delta <- eval[k]
    for (i in 1:d_size){
      dl <- D_l[i]
      d_u <- dl / (delta * dl + 1)
      d_e <- d_u * delta
      OmegaU[i, k] <- d_u
      OmegaE[i, k] <- d_e
    }
  }
  return(list(OmegaU, OmegaE))
}
