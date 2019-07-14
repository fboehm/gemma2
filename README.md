
## Statuses

[![Travis-CI Build
Status](https://travis-ci.org/fboehm/gemma2.svg?branch=master)](https://travis-ci.org/fboehm/gemma2)

[![codecov](https://codecov.io/gh/fboehm/gemma2/branch/master/graph/badge.svg)](https://codecov.io/gh/fboehm/gemma2)

<!-- README.md is generated from README.Rmd. Please edit that file -->

## Overview

`gemma2` is an implementation in R of the GEMMA v 0.97 EM algorithm that
is part of the GEMMA algorithm for REML estimation of multivariate
linear mixed effects models of the form:

\[vec(Y) = X vec(B) + vec(G) + vec(E)\]

where \(E\) is a n by 2 matrix of random effects that follows the
matrix-variate normal distribution

\[G \sim MN(0, K, V_g)\]

where \(K\) is a relatedness matrix and \(V_g\) is a 2 by 2 covariance
matrix for the two traits of interest.

Additionally, the random errors matrix \(E\) follows the distribution:

\[E \sim MN(0, I_n, V_e)\]

and \(G\) and \(E\) are independent.

## Installation

To install `gemma2`, use the `devtools` R package from CRAN. If you
haven’t installed `devtools`, please run this line of code:

``` r
install.packages("devtools")
```

Then, run this line of code to install `gemma2`:

``` r
devtools::install_github("fboehm/gemma2")
```

## References

X. Zhou & M. Stephens. Efficient multivariate linear mixed model
algorithms for genome-wide association studies. *Nature Methods* volume
11, pages 407–409 (2014). <https://www.nature.com/articles/nmeth.2848>

<https://github.com/genetics-statistics/GEMMA>
