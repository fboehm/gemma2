#' Calculate log likelihood
#'
#' @param eval eigenvalues vector from decomposition of relatedness matrix
#' @param D_l vector
#' @param Qi inverse of Q matrix
#' @param UltVehiY matrix of (transformed) Y values
#' @param xHiy vector
#' @export
MphCalcLogL <- function(eval, D_l, Qi, UltVehiY, xHiy){
  n_size <- length(eval)
  d_size <- length(D_l)
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


#' Perform EM algorithm
#'
#' @param func_name indicator for calculating restricted or unrestricted likelihood
#' @param max_iter maximum number of iterations for EM algorithm
#' @param max_prec maximum precision for EM algorithm
#' @param eval vector of eigenvalues from relatedness matrix decomposition
#' @param X design matrix
#' @param Y matrix of phenotype values
#' @param V_g covariance matrix
#' @param V_e covariance matrix
#' @export
MphEM <- function(func_name = "R", max_iter = 10000, max_prec = 1 / 1000,
                  eval, X, Y, V_g, V_e){
  n_size <- length(eval)
  c_size <- nrow(X)
  d_size <- nrow(Y)
  dc_size <- c_size * d_size
  # calculate XXt and XXti
  XXt <- X %*% t(X)
  #XXti <- solve(XXt)
  ## XXti is only used for unrestricted likelihood!
  # calculate constant for logl
  #if (func_name == "R"){
  logl_const <- (n_size - c_size) * d_size * log(2 * pi) / 2 + d_size * log(det(XXt)) / 2
  print(logl_const)
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
    #if (t > 1 & abs(logl_new - logl_old) < max_prec) {break}
    logl_old <- logl_new
    co_out <- calc_omega(eval, D_l)
    co_out[[1]] -> OmegaU
    co_out[[2]] -> OmegaE
    if (func_name == 'R') {
      UltVehiB <- UpdateRL_B(xHiy, Qi)
      UltVehiBX <- UltVehiB %*% X
    }
    UltVehiU <- update_u(OmegaE, UltVehiY, UltVehiBX)
    UltVehiE <- update_e(UltVehiY, UltVehiBX, UltVehiU)
    U_hat <- t(UltVeh) %*% UltVehiU
    E_hat <- t(UltVeh) %*% UltVehiE
    B <- t(UltVeh) %*% UltVehiB
    cs_out <- calc_sigma(func_name = func_name, eval = eval, D_l = D_l, 
                         X = X, OmegaU = OmegaU, OmegaE = OmegaE, 
                         UltVeh = UltVeh, Qi = Qi)
    cs_out[[1]] -> Sigma_uu
    cs_out[[2]] -> Sigma_ee
    uv_out <- update_v(eval, U_hat, E_hat, Sigma_uu, Sigma_ee)
    # update V_g and V_e
    uv_out[[1]] -> V_g
    uv_out[[2]] -> V_e
    out[[t]] <- list(logl_new, V_g, V_e, Sigma_uu, Sigma_ee, B, 
                     U_hat, E_hat, OmegaU, OmegaE, logdet_Ve, UltVeh, UltVehi, 
                     D_l, xHiy)
  }
  return(out)
}
