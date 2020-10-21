#' Stagger matrices within a larger, block-diagonal matrix
#'
#' @param ... one or more matrices, separated by commas
#' @return a block-diagonal matrix, with the inputted matrices as blocks on the diagonal.
#' @examples
#' foo <- matrix(rnorm(40000), ncol = 8)
#' block_diag <- stagger_mats(foo, foo)
#' dim(foo)
#' dim(block_diag)
#' @export

stagger_mats <- function(...){
  return(bdiag_m(...))
}


bdiag_m <- function(lmat) {
  ## Copyright (C) 2016 Martin Maechler, ETH Zurich
  if(!length(lmat)) return(new("dgCMatrix"))
  stopifnot(is.list(lmat), is.matrix(lmat[[1]]),
            (k <- (d <- dim(lmat[[1]]))[1]) == d[2], # k x k
            all(vapply(lmat, dim, integer(2)) == k)) # all of them
  N <- length(lmat)
  if(N * k > .Machine$integer.max)
    stop("resulting matrix too large; would be  M x M, with M=", N*k)
  M <- as.integer(N * k)
  ## result: an   M x M  matrix
  new("dgCMatrix", Dim = c(M,M),
      ## 'i :' maybe there's a faster way (w/o matrix indexing), but elegant?
      i = as.vector(matrix(0L:(M-1L), nrow=k)[, rep(seq_len(N), each=k)]),
      p = k * 0L:M,
      x = as.double(unlist(lmat, recursive=FALSE, use.names=FALSE)))
}
