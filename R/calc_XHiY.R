#' Calculate XHiY
#'
#' @param eval vector of eigenvalues from the decomposition of the relatedness matrix
#' @param D_l vector of length d_size
#' @param X design matrix
#' @param UltVehiY a matrix
#' @return numeric vector
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
#' eigen2(kinship) -> eout
#' eout$values -> eval
#' eout$vectors -> U
#' UltVehi <- matrix(c(0, -1.76769, -1.334414, 0),
#' nrow = 2,
#' byrow = FALSE) # from output of eigen_proc()
#' calc_XHiY(eval = eval,
#' D_l = c(0.9452233, 5.9792268),
#'           X = rep(1, 100) %*% U,
#'           UltVehiY = UltVehi %*% t(phe16) %*% U
#'           )


calc_XHiY <- function(eval, D_l, X, UltVehiY){
  # check inputs
  stopifnot(length(eval) == ncol(X),
            is.vector(eval),
            is.vector(D_l))
  n_size <- length(eval)
  c_size <- nrow(X)
  d_size <- length(D_l)
  xHiy <- rep(0, d_size * c_size)
  for (i in 1:d_size){
    dl <- D_l[i]
    for (j in 1:c_size){
      d <- 0
      for (k in 1:n_size){
        x <- X[j, k]
        y <- UltVehiY[i, k]
        delta <- eval[k]
        d <- d + x * y / (delta * dl + 1)
      }
      xHiy[(j - 1) * d_size + i] <- d
    }
  }
  return(xHiy)
}


#' Eigendecomposition procedure for Vg and Ve
#'
#' @param V_g a d_size by d_size covariance matrix
#' @param V_e a d_size by d_size covariance matrix
#' @param tol a positive number indicating the tolerance for isSymmetric
eigen_proc <- function(V_g, V_e, tol = 1 / 10000){
  # check inputs
  stopifnot(isSymmetric(V_g, tol = tol), isSymmetric(V_e, tol = tol))
  d_size <- nrow(V_g)
  logdet_Ve <- 0
  Lambda <- matrix(nrow = d_size, ncol = d_size)
  V_e_temp <- matrix(nrow = d_size, ncol = d_size)
  V_e_h <- matrix(0, nrow = d_size, ncol = d_size)
  V_e_hi <- matrix(0, nrow = d_size, ncol = d_size)
  VgVehi <- matrix(nrow = d_size, ncol = d_size)
  U_l <- matrix(nrow = d_size, ncol = d_size)
  V_e -> V_e_temp
  eigen2(V_e_temp) -> eout
  eout$values -> D_l
  eout$vectors -> U_l
  if (length(U_l == 1)) U_l <- as.matrix(U_l)
  for (i in 1:d_size){
    d <- D_l[i]
    if (d > 0){
      logdet_Ve <- logdet_Ve + log(d)
      U_col <- U_l[, i]
      d <- sqrt(d)
      V_e_h <- V_e_h + d * U_col %*% t(U_col)
      V_e_hi <- V_e_hi + U_col %*% t(U_col) / d
    }
  }
  V_g %*% V_e_hi -> VgVehi
  Lambda <- V_e_hi %*% VgVehi

  eigen2(Lambda) -> eout
  eout$values -> D_l
  eout$vectors -> U_l
  if (length(U_l) == 1) U_l <- as.matrix(U_l)
  D_l[D_l < 0] <- 0
  UltVeh <- t(U_l) %*% V_e_h
  UltVehi <- t(U_l) %*% V_e_hi

  return(list(logdet_Ve = logdet_Ve, UltVeh = UltVeh, UltVehi = UltVehi, D_l = D_l ))
}

#' Calculate Qi (inverse of Q) and log determinant of Q
#'
#' @param eval vector of eigenvalues from decomposition of relatedness matrix
#' @param D_l vector of length d_size
#' @param X design matrix
#' @return a list of length two. First entry in the list is a symmetric numeric matrix, Qi, the inverse of the Q matrix. The second entry in the outputted list is the log determinant of the matrix Q for use in likelihood calculations.
#' @examples
#' as.matrix(readr::read_tsv(system.file("extdata",
#' "mouse100.cXX.txt",
#' package = "gemma2"),
#' col_names = FALSE)[, 1:100]) -> kinship
#' eigen2(kinship) -> e2_out
#' e2_out$values -> eval
#' e2_out$vectors -> U
#' eigen_proc(V_g = diag(c(1.91352, 0.530827)),
#' V_e = diag(c(0.320028, 0.561589))) -> ep_out
#'
#' calc_qi(eval = eval,
#' D_l = ep_out[[4]],
#' X = t(rep(1, 100)) %*% U)
calc_qi <- function(eval, D_l, X){
  n_size <- length(eval)
  d_size <- length(D_l)
  c_size <- nrow(X)# what is c_size? it's the number of rows in the transposed genotypes matrix
  # transposed genotypes matrix is c by n
  dc_size <- d_size * c_size
  Q <- matrix(0, nrow = dc_size, ncol = dc_size)
  # What is this loop actually calculating?
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
#' @return list of length 2. First entry in the list is the symmetric matrix OmegaU. Second entry in the list is the symmetric matrix OmegaE.
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
