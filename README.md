
-   [Statuses](#statuses)
-   [Overview](#overview)
-   [Installation](#installation)

Statuses
--------

[![Travis-CI Build Status](https://travis-ci.org/fboehm/gemma2.svg?branch=master)](https://travis-ci.org/fboehm/gemma2)

[![codecov](https://codecov.io/gh/fboehm/gemma2/branch/master/graph/badge.svg)](https://codecov.io/gh/fboehm/gemma2)

<!-- README.md is generated from README.Rmd. Please edit that file -->
Overview
--------

`gemma2` is an implementation in R of the GEMMAv 0.97 EM algorithm that is part of the GEMMA algorithm for REML estimation of multivariate linear mixed effects models of the form:

*v**e**c*(*Y*)=*X**v**e**c*(*B*)+*v**e**c*(*G*)+*v**e**c*(*E*)

where *E* is a n by 2 matrix of random effects that follows the matrix-variate normal distribution

*G* ∼ *M**N*(0, *K*, *V*<sub>*g*</sub>)

where *K* is a relatedness matrix and *V*<sub>*g*</sub> is a 2 by 2 covariance matrix for the two traits of interest.

Additionally, the random errors matrix *E* follows the distribution:

*E* ∼ *M**N*(0, *I*<sub>*n*</sub>, *V*<sub>*e*</sub>)

and *G* and *E* are independent.

Installation
------------

To install `gemma2`, use the `devtools` R package from CRAN. If you haven't installed `devtools`, please run this line fo code:

``` r
install.packages("devtools")
```

Then, run this line of code to install `gemma2`:

``` r
devtools::install_github("fboehm/gemma2")
```
