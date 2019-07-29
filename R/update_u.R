#' Update U matrix
#'
#' @param OmegaE the OmegaE matrix, calculated in calc_omega
#' @param UltVehiY matrix
#' @param UltVehiBX matrix
#' @export
#' @family expectation-maximization functions
#' @examples
#' readr::read_tsv(system.file("extdata", "mouse100.pheno.txt", package = "gemma2"), col_names = FALSE) -> pheno
#' phe16 <- as.matrix(pheno[, c(1, 6)])
#' as.matrix(readr::read_tsv(system.file("extdata", "mouse100.cXX.txt", package = "gemma2"), col_names = FALSE)[, 1:100]) -> kinship
#' eigen2(kinship) -> e2_out
#' e2_out$values -> eval
#' e2_out$vectors -> U
#' eigen_proc(V_g = diag(c(1.91352, 0.530827)), V_e = diag(c(0.320028, 0.561589))) -> ep_out
#' UltVehi <- ep_out[[3]]
#' calc_omega(eval, ep_out$D_l) -> co_out
#' update_u(OmegaE = co_out[[2]],
#'         UltVehiY = UltVehi %*% t(phe16),
#'         UltVehiBX = matrix(c(-0.71342, -0.824482), ncol = 1) %*% t(rep(1, 100))
#' )
update_u <- function(OmegaE, UltVehiY, UltVehiBX){
  #void UpdateU (const gsl_matrix *OmegaE, const gsl_matrix *UltVehiY, const gsl_matrix *UltVehiBX, gsl_matrix *UltVehiU)
  #{
  #  gsl_matrix_memcpy (UltVehiU, UltVehiY);
  UltVehiU <- UltVehiY
  #  gsl_matrix_sub (UltVehiU, UltVehiBX);
  UltVehiU <- UltVehiU - UltVehiBX
  #
  #  gsl_matrix_mul_elements (UltVehiU, OmegaE);
  UltVehiU <- UltVehiU * OmegaE
  #  return;
  #}
  #UltVehiU <- (UltVehiY - UltVehiBX) * OmegaE
  return(UltVehiU)
}

#' Update E
#'
#' @param UltVehiY matrix of transformed Y values
#' @param UltVehiBX matrix of transformed BX values
#' @param UltVehiU matrix of transformed U values
#' @family expectation-maximization functions
#' @export
update_e <- function(UltVehiY, UltVehiBX, UltVehiU){
  #void UpdateE (const gsl_matrix *UltVehiY, const gsl_matrix *UltVehiBX, const gsl_matrix *UltVehiU, gsl_matrix *UltVehiE)
  #{
   # gsl_matrix_memcpy (UltVehiE, UltVehiY);
  #UltVehiY -> UltVehiE
  #  gsl_matrix_sub (UltVehiE, UltVehiBX);
  #UltVehiE <- UltVehiE - UltVehiBX
  #  gsl_matrix_sub (UltVehiE, UltVehiU);
  #UltVehiE <- UltVehiE - UltVehiU
  #
  #  return;
  #}
  UltVehiE <- UltVehiY - UltVehiBX - UltVehiU
  return(UltVehiE)
}

#' Update B for restricted log likelihood
#'
#' @param xHiy vector
#' @param Qi Q inverse matrix
#' @param d_size number of traits
#' @family expectation-maximization functions
#' @export
UpdateRL_B <- function(xHiy, Qi, d_size){
  #void UpdateRL_B (const gsl_vector *xHiy, const gsl_matrix *Qi, gsl_matrix *UltVehiB)
  #{
  #  size_t d_size=UltVehiB->size1, c_size=UltVehiB->size2, dc_size=Qi->size1;
  nrow(Qi) -> dc_size
  c_size <- dc_size / d_size
  #
  #  gsl_vector *b=gsl_vector_alloc (dc_size);
  #
  #  //calculate b=Qiv
  #  gsl_blas_dgemv(CblasNoTrans, 1.0, Qi, xHiy, 0.0, b);
  b <- Qi %*% xHiy
  UltVehiB <- matrix(nrow = d_size, ncol = c_size)
  #
  #  //copy b to UltVehiB
  #  for (size_t i=0; i<c_size; i++) {
  for (i in 1:c_size){
    #    gsl_vector_view UltVehiB_col=gsl_matrix_column (UltVehiB, i);
    b_subcol <- b[(1 + (i - 1) * d_size):(i * d_size)]
  #    gsl_vector_const_view b_subcol=gsl_vector_const_subvector (b, i*d_size, d_size);
  #    gsl_vector_memcpy (&UltVehiB_col.vector, &b_subcol.vector);
  #  }
    b_subcol -> UltVehiB[, i] # could use as.matrix here
  }#
  #  gsl_vector_free(b);
  #
  #  return;
  #}
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
#' @export
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
#' @export
calc_sigma <- function(eval, D_l, X, OmegaU, OmegaE, UltVeh, Qi){
  #void CalcSigma (const char func_name, const gsl_vector *eval, const gsl_vector *D_l, const gsl_matrix *X, const gsl_matrix *OmegaU, const gsl_matrix *OmegaE, const gsl_matrix *UltVeh, const gsl_matrix *Qi, gsl_matrix *Sigma_uu, gsl_matrix *Sigma_ee)
  #{
  #  if (func_name!='R' && func_name!='L' && func_name!='r' && func_name!='l') {cout<<"func_name only takes 'R' or 'L': 'R' for log-restricted likelihood, 'L' for log-likelihood."<<endl; return;}
  #
  #  size_t n_size=eval->size, c_size=X->size1, d_size=D_l->size, dc_size=Qi->size1;
  n_size <- length(eval)
  c_size <- nrow(X)
  d_size <- length(D_l)
  dc_size <- nrow(Qi)
  #  gsl_matrix_set_zero(Sigma_uu);
  #  gsl_matrix_set_zero(Sigma_ee);
  Sigma_ee <- matrix(0, nrow = d_size, ncol = d_size)
  Sigma_uu <- Sigma_ee
  #
  #  double delta, dl, x, d;
  #
   # //calculate the first diagonal term
  #  gsl_vector_view Suu_diag=gsl_matrix_diagonal (Sigma_uu);
  #  gsl_vector_view See_diag=gsl_matrix_diagonal (Sigma_ee);
  #
  #  for (size_t k=0; k<n_size; k++) {
  for (k in 1:n_size){
    #    gsl_vector_const_view OmegaU_col=gsl_matrix_const_column (OmegaU, k);
    OmegaU_col <- OmegaU[, k]
  #    gsl_vector_const_view OmegaE_col=gsl_matrix_const_column (OmegaE, k);
    OmegaE_col <- OmegaE[, k]
  #
  #    gsl_vector_add (&Suu_diag.vector, &OmegaU_col.vector);

    diag(Sigma_uu) <- diag(Sigma_uu) + OmegaU_col
  #    gsl_vector_add (&See_diag.vector, &OmegaE_col.vector);
    diag(Sigma_ee) <- diag(Sigma_ee) + OmegaE_col
  }
  #  }
  #
  #  //calculate the second term for reml
  #  if (func_name=='R' || func_name=='r') {
  #    gsl_matrix *M_u=gsl_matrix_alloc(dc_size, d_size);
  #    gsl_matrix *M_e=gsl_matrix_alloc(dc_size, d_size);
  #    gsl_matrix *QiM=gsl_matrix_alloc(dc_size, d_size);
  #
  #    gsl_matrix_set_zero(M_u);
  M_u <- matrix(0, nrow = dc_size, ncol = d_size)
  M_e <- M_u
  #    gsl_matrix_set_zero(M_e);
  #
  #    for (size_t k=0; k<n_size; k++) {
  #      delta=gsl_vector_get(eval, k);
  #      //if (delta==0) {continue;}
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
    #
  #      for (size_t i=0; i<d_size; i++) {
  #        dl=gsl_vector_get(D_l, i);
  #        for (size_t j=0; j<c_size; j++) {
  #          x=gsl_matrix_get(X, j, k);
  #          d=x/(delta*dl+1.0);
  #          gsl_matrix_set(M_e, j*d_size+i, i, d);
  #          gsl_matrix_set(M_u, j*d_size+i, i, d*dl);
  #        }
  #      }
  #      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Qi, M_u, 0.0, QiM);
    QiM <- Qi %*% M_u
      #      gsl_blas_dgemm(CblasTrans, CblasNoTrans, delta, M_u, QiM, 1.0, Sigma_uu);
    Sigma_uu <- Sigma_uu + t(M_u) %*% QiM * delta
    QiM <- Qi %*% M_e
    Sigma_ee <- Sigma_ee + t(M_e) %*% QiM
  } # end loop in k
    #
  #      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Qi, M_e, 0.0, QiM);
  #      gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, M_e, QiM, 1.0, Sigma_ee);
  #    }
  #
  #    gsl_matrix_free(M_u);
  #    gsl_matrix_free(M_e);
  #    gsl_matrix_free(QiM);
  # }
  #
  #  //multiply both sides by VehUl
  #  gsl_matrix *M=gsl_matrix_alloc (d_size, d_size);
  #
  #  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Sigma_uu, UltVeh, 0.0, M);
  #  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, UltVeh, M, 0.0, Sigma_uu);
  #  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Sigma_ee, UltVeh, 0.0, M);
  #  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, UltVeh, M, 0.0, Sigma_ee);
  #
  #  gsl_matrix_free(M);
  #  return;
  #}
  M <- Sigma_uu %*% UltVeh
  Sigma_uu <- t(UltVeh) %*% M
  M <- Sigma_ee %*% UltVeh
  Sigma_ee <- t(UltVeh) %*% M
  return(list(Sigma_ee, Sigma_uu))
}

