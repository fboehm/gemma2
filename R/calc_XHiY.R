#' Calculate XHiY
#'
#' @param eval vector of eigenvalues from the decomposition of the relatedness matrix
#' @param D_l vector of length d_size
#' @param X design matrix
#' @param UltVehiY a matrix
#' @export
calc_XHiY <- function(eval, D_l, X, UltVehiY){
  # check inputs
  stopifnot(length(eval) == ncol(X),
            is.vector(eval),
            is.vector(D_l))
  n_size <- length(eval)
  c_size <- nrow(X) - 1
  d_size <- length(D_l)
  xHiy <- numeric()
  for (i in 1:d_size){
    dl <- D_l[i]
    for (j in (1 + 1):(c_size + 1)){
      d <- 0
      for (k in 1:n_size){
        x <- X[j, k]
        y <- UltVehiY[i, k]
        delta <- eval[k]
        d <- d + x * y / (delta * dl + 1)
      }
      xHiy[(j - 1 - 1) * d_size + i] <- d
    }
  }
  return(xHiy)
}

#' Eigendecomposition procedure for Vg and Ve
#'
#' @param V_g a d_size by d_size covariance matrix
#' @param V_e a d_size by d_size covariance matrix
#' @export
#double EigenProc (const gsl_matrix *V_g, const gsl_matrix *V_e, gsl_vector *D_l, gsl_matrix *UltVeh, gsl_matrix *UltVehi)
##{
eigen_proc <- function(V_g, V_e){
  # check inputs
  stopifnot(isSymmetric(V_g), isSymmetric(V_e))
#  size_t d_size=V_g->size1;
  d_size <- nrow(V_g)
#  double d, logdet_Ve=0.0;
  logdet_Ve <- 0
#
#  //eigen decomposition of V_e
#  gsl_matrix *Lambda=gsl_matrix_alloc (d_size, d_size);
  Lambda <- matrix(nrow = d_size, ncol = d_size)
#  gsl_matrix *V_e_temp=gsl_matrix_alloc (d_size, d_size);
  V_e_temp <- matrix(nrow = d_size, ncol = d_size)
#  gsl_matrix *V_e_h=gsl_matrix_alloc (d_size, d_size);
  V_e_h <- matrix(0, nrow = d_size, ncol = d_size)
#  gsl_matrix *V_e_hi=gsl_matrix_alloc (d_size, d_size);
  V_e_hi <- matrix(0, nrow = d_size, ncol = d_size)
#  gsl_matrix *VgVehi=gsl_matrix_alloc (d_size, d_size);
  VgVehi <- matrix(nrow = d_size, ncol = d_size)
#  gsl_matrix *U_l=gsl_matrix_alloc (d_size, d_size);
  U_l <- matrix(nrow = d_size, ncol = d_size)
#
#  gsl_matrix_memcpy(V_e_temp, V_e);
  V_e -> V_e_temp
#  EigenDecomp(V_e_temp, U_l, D_l, 0);
  eigen(V_e_temp) -> eout
  eout$values -> D_l
  eout$vectors -> U_l
#
#  //calculate V_e_h and V_e_hi
#  gsl_matrix_set_zero(V_e_h);
#  gsl_matrix_set_zero(V_e_hi);
#  for (size_t i=0; i<d_size; i++) {
  for (i in 1:d_size){
#    d=gsl_vector_get (D_l, i);
    d <- D_l[i]
#    if (d<=0) {continue;}
    if (d > 0){
#    logdet_Ve+=log(d);
      logdet_Ve <- logdet_Ve + log(d)
#
#    gsl_vector_view U_col=gsl_matrix_column(U_l, i);
      U_col <- U_l[, i]
      d <- sqrt(d)
#    d=sqrt(d);
#    gsl_blas_dsyr (CblasUpper, d, &U_col.vector, V_e_h);
      V_e_h <- V_e_h + d * U_col %*% t(U_col)
#    d=1.0/d;
#    gsl_blas_dsyr (CblasUpper, d, &U_col.vector, V_e_hi);
      V_e_hi <- V_e_hi + U_col %*% t(U_col) / d
    }
  }
#  }
#
#  //copy the upper part to lower part
#  for (size_t i=0; i<d_size; i++) {
#    for (size_t j=0; j<i; j++) {
#      gsl_matrix_set (V_e_h, i, j, gsl_matrix_get(V_e_h, j, i));
#      gsl_matrix_set (V_e_hi, i, j, gsl_matrix_get(V_e_hi, j, i));
#    }
#  }
#
#  //calculate Lambda=V_ehi V_g V_ehi
#  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, V_g, V_e_hi, 0.0, VgVehi);
  V_g %*% V_e_hi -> VgVehi

#  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, V_e_hi, VgVehi, 0.0, Lambda);
  Lambda <- V_e_hi %*% VgVehi

    #
#  //eigen decomposition of Lambda
#  EigenDecomp(Lambda, U_l, D_l, 0);
  eigen(Lambda) -> eout
  eout$values -> D_l
  eout$vectors -> U_l

    #
#  for (size_t i=0; i<d_size; i++) {
#    d=gsl_vector_get (D_l, i);
  D_l[D_l < 0] <- 0
#    if (d<0) {gsl_vector_set (D_l, i, 0);}
#  }
#
#  //calculate UltVeh and UltVehi
#  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, U_l, V_e_h, 0.0, UltVeh);
  UltVeh <- t(U_l) %*% V_e_h
#  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, U_l, V_e_hi, 0.0, UltVehi);
  UltVehi <- t(U_l) %*% V_e_hi

#
#    //free memory
#  gsl_matrix_free (Lambda);
#  gsl_matrix_free (V_e_temp);
#  gsl_matrix_free (V_e_h);
#  gsl_matrix_free (V_e_hi);
#  gsl_matrix_free (VgVehi);
#  gsl_matrix_free (U_l);
#
  return(list(logdet_Ve, UltVeh, UltVehi, D_l))
}
#  return logdet_Ve;
#}
#' Calculate Qi (inverse of Q) and log determinant of Q
#'
#' @param eval vector of eigenvalues from decomposition of relatedness matrix
#' @param D_l vector of length d_size
#' @param X a design matrix
#' @export
calc_qi <- function(eval, D_l, X){
  n_size <- length(eval)
  d_size <- length(D_l)
  c_size <- nrow(X) - 1# what is c_size? c + 1 is the row number for the X design matrix
  dc_size <- d_size * c_size
  Q <- matrix(0, nrow = dc_size, ncol = dc_size)

  for (i in (1 + 1):(c_size +1)){
    for (j in (1 + 1):(c_size + 1)){
      for (l in 1:d_size){
        dl <- D_l[l]
        if (j < i){
          d <- Q[(j - 1 - 1) * d_size + l, (i - 1 - 1) * d_size + l]
        } else {
          d <- 0
          for (k in 1:n_size){
            d1 <- X[i, k]
            d2 <- X[j, k]
            delta <- eval[k]
            d <- d + d1 * d2 / (dl * delta + 1)
          }
        }
        Q[(i - 1 - 1) * d_size + l, (j - 1 - 1) * d_size + l] <- d
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
