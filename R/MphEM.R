#' Perform EM algorithm
#'
#' @param max_iter maximum number of iterations for EM algorithm
#' @param max_prec maximum precision for EM algorithm
#' @param eval vector of eigenvalues from relatedness matrix decomposition
#' @param X design matrix
#' @param Y matrix of phenotype values
#' @param V_g genetic covariance matrix
#' @param V_e error covariance matrix
#' @return a list of lists. Length of list corresponds to number of EM iterations
#' @export
MphEM <- function(max_iter = 10000, max_prec = 1 / 10000,
                  eval, X, Y, V_g, V_e){
  n_size <- length(eval)
  # be careful with defining c_size here

  c_size <- nrow(X)
  d_size <- nrow(Y)
  dc_size <- c_size * d_size
  # calculate XXt and XXti
  XXt <- X %*% t(X)
  #XXti <- solve(XXt)
  ## XXti is only used for unrestricted likelihood!
  # calculate constant for logl
  #if (func_name == "R"){
  logl_const <- - (n_size - c_size) * d_size * log(2 * pi) / 2 + d_size * log(det(XXt)) / 2
  #}
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
    cq_out <- calc_qi(eval, D_l, X)
    cq_out[[1]] -> Qi
    cq_out[[2]] -> lndetQ
    #
    UltVehiY <- UltVehi %*% Y
    #
    xHiy <- calc_XHiY(eval, D_l, X, UltVehiY)
    logl_new <- logl_const + MphCalcLogL(eval = eval, xHiy = xHiy,
                                         D_l = D_l, UltVehiY = UltVehiY,
                                         Qi = Qi) - 0.5 * n_size * logdet_Ve
    #if (func_name=='R') {
    logl_new <- logl_new - 0.5 * (lndetQ - c_size * logdet_Ve)
    #}
    if (t > 1){
      if (logl_new - logl_old < max_prec){break}
    }
    #if (t > 1 & abs(logl_new - logl_old) < max_prec) {break}
    logl_old <- logl_new
    co_out <- calc_omega(eval, D_l)
    co_out[[1]] -> OmegaU
    co_out[[2]] -> OmegaE
    UltVehiB <- UpdateRL_B(xHiy, Qi, d_size = nrow(Y))
    UltVehiBX <- UltVehiB %*% X

    UltVehiU <- update_u(OmegaE, UltVehiY, UltVehiBX)
    UltVehiE <- update_e(UltVehiY, UltVehiBX, UltVehiU)
    U_hat <- t(UltVeh) %*% UltVehiU
    E_hat <- t(UltVeh) %*% UltVehiE
    B <- t(UltVeh) %*% UltVehiB
    cs_out <- calc_sigma(eval = eval, D_l = D_l,
                         X = X, OmegaU = OmegaU, OmegaE = OmegaE,
                         UltVeh = UltVeh, Qi = Qi)
    cs_out[[1]] -> Sigma_ee
    cs_out[[2]] -> Sigma_uu
    uv_out <- update_v(eval, U_hat, E_hat, Sigma_uu, Sigma_ee)
    # update V_g and V_e
    uv_out[[1]] -> V_g
    uv_out[[2]] -> V_e
    out[[t]] <- list(logl_new = logl_new, Vg = V_g, Ve = V_e, Sigma_uu = Sigma_uu,
                     Sigma_ee = Sigma_ee, B = B,
                     U_hat = U_hat, E_hat = E_hat,
                     OmegaU = OmegaU, OmegaE = OmegaE, logdet_Ve = logdet_Ve,
                     UltVeh = UltVeh, UltVehi = UltVehi,
                     Dl = D_l, xHiy = xHiy, logl_const = logl_const, UltVehiU = UltVehiU
                     )
  }
  return(out)
}

