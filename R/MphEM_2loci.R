#' Perform EM algorithm for pleiotropy or close linkage
#'
#' @param max_iter maximum number of iterations for EM algorithm
#' @param max_prec maximum precision for EM algorithm
#' @param eval vector of eigenvalues from relatedness matrix decomposition
#' @param X1 transformed genotype probabilities matrix 1
#' @param X2 transformed genotype probabilities matrix 2
#' @param Y transformed matrix of phenotype values
#' @param V_g genetic covariance matrix
#' @param V_e error covariance matrix
#' @export
MphEM_2loci <- function(max_iter = 10000, max_prec = 1 / 1000, eval, X1, X2, Y, V_g = diag(1, 2), V_e = diag(1, 2)){
  n_size <- length(eval)
  # be careful with defining c_size here
  stopifnot(nrow(X1) == nrow(X2), ncol(X1) == length(eval))
  c_size <- nrow(X1)
  d_size <- nrow(Y)
  dc_size <- c_size * d_size
  # convert Y to vector Yv
  Yv <- t(as.vector(t(Y)))
  # stagger X1 & X2
  X <- t(pleiotropy::stagger_mats(t(X1), t(X2)))
  # calculate XXt and XXti
  XXt <- X %*% t(X)
  #XXti <- solve(XXt)
  ## XXti is only used for unrestricted likelihood!
  # calculate constant for logl for REML
  logl_const <- - (n_size - c_size) * d_size * log(2 * pi) / 2 + d_size * log(sqrt(det(XXt))) / 2
  # NOTE USE OF SQRT ABOVE!
  out <- list()
  #logl_old <- 1 # we need to define logl_old before using it within the EM iterations
  # start EM
  for (t in 1:max_iter){
    # eigen_proc
    ep_out <- eigen_proc(V_g, V_e)
    ep_out[[1]] -> logdet_Ve
    ep_out[[2]] -> UltVeh
    ep_out[[3]] -> UltVehi
    ep_out[[4]] -> D_l
    # calc_qi
    cq_out <- calc_qi_2loci(eval, D_l, X1, X2)
    cq_out[[1]] -> Qi
    cq_out[[2]] -> lndetQ
    #
    UltVehiY <- UltVehi %*% Y
    #
    xHiy <- calc_XHiY_2loci(eval, D_l, X1, X2, UltVehiY)

    logl_new <- logl_const + MphCalcLogL(eval = eval, xHiy = xHiy,
                                         D_l = D_l, UltVehiY = UltVehiY,
                                         Qi = Qi) - 0.5 * n_size * logdet_Ve
    #if (func_name=='R') {
    logl_new <- logl_new - 0.5 * (lndetQ - c_size * logdet_Ve)#}
    if (t > 1){
      if (logl_new - logl_old < max_prec){break}
    }
    #if (t > 1 & abs(logl_new - logl_old) < max_prec) {break}
    logl_old <- logl_new
    co_out <- calc_omega(eval, D_l)
    co_out[[1]] -> OmegaU
    co_out[[2]] -> OmegaE
    UltVehiB <- UpdateRL_B(xHiy, Qi, d_size = nrow(Y))
    #UltVehiBX <- UltVehiB %*% X
    #foo <- UltVehiB %*% X # This gives repeated entries.

    #foo <- t(X) %*% as.vector(UltVehiB)
    #UltVehiBX <- rbind(foo[ 1:n_size, ], foo[(n_size + 1):(2 * n_size), ])
    t(as.matrix(as.vector(t(UltVehiB)))) -> UltVehiB_vec
    UltVehiB_vec %*% X -> UltVehiB_vecX
    ## fix dimensionality of UltVehiBX to be 2 rows by n_size columns, where
    ## n_size is the number of mice
    UltVehiBX <- matrix(nrow = 2, ncol = n_size, data = as.vector(UltVehiB_vecX), byrow = TRUE)
    # Check this (above) !! DONE.

    # apply update_u to both UltVehiBX1 and UltVehiBX2
    UltVehiU <- update_u(OmegaE, UltVehiY, UltVehiBX)
    # Now, look at update_e()
    UltVehiE <- update_e(UltVehiY, UltVehiBX, UltVehiU)
    #
    U_hat <- t(UltVeh) %*% UltVehiU


    E_hat <- t(UltVeh) %*% UltVehiE
    B <- t(UltVeh) %*% UltVehiB
    cs_out <- calc_sigma_2loci(eval = eval, D_l = D_l,
                               X1 = X1, X2 = X2, OmegaU = OmegaU, OmegaE = OmegaE,
                               UltVeh = UltVeh, Qi = Qi)
    cs_out[[1]] -> Sigma_ee
    cs_out[[2]] -> Sigma_uu
    uv_out <- update_v(eval, U_hat, E_hat, Sigma_uu, Sigma_ee)
    # update V_g and V_e
    uv_out[[1]] -> V_g
    uv_out[[2]] -> V_e
    out[[t]] <- list(logl_new, V_g, V_e, Sigma_uu, Sigma_ee, B,
                     U_hat, E_hat, OmegaU, OmegaE, logdet_Ve, UltVeh, UltVehi,
                     D_l, xHiy, logl_const, UltVehiU)
  }
  return(out)
}


