## Test environments
* local OS X install, R 3.6.1
* ubuntu 14.04 (on travis-ci), R 3.5.3
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 notes

This is a resubmission. I previously submitted gemma2 v 0.0.2. 

## Replies to CRAN submission feedback

On 14.07.2019 19:55, CRAN submission wrote:
> [This was generated from CRAN.R-project.org/submit.html]
> 
> The following package was uploaded to CRAN:
> ===========================================
> 
> Package Information:
> Package: gemma2
> Version: 0.0.2
> Title: Zhou & Stephens (2014) GEMMA Multivariate Linear Mixed Model


### Point 1

> Thanks, please write references in your Description text and not in your 
> title. Maybe just:
> GEMMA Multivariate Linear Mixed Model

Thank you for the recommendation. I've made this change to the title.

> Author(s): Frederick Boehm [aut, cre]
>    (<https://orcid.org/0000-0002-1644-5931>)
> Maintainer: Frederick Boehm <frederick.boehm@gmail.com>
> Suggests: covr, testthat, knitr, rmarkdown
> Description: Fits a multivariate linear mixed effects model that uses a
>    polygenic term, after Zhou & Stephens (2014)
>    (<https://www.nature.com/articles/nmeth.2848>). Of particular
>    interest is the estimation of variance components with REML.


### Point 2

> Please explain all acronyms (e.g. GEMMA) in your Description text to 
avoid misunderstandings.

Thanks. I've added text in the Description field to explain GEMMA. I've avoided using acronyms in the Description field.


### Point 3

> Please add more examples in your Rd-files. Since you have tests, you can 
> wrap them in \donttest{}.

Thank you. I've added examples for all exported functions.

### Point 4

> Missing Rd-tags:
  calc_omega.Rd: \value
  calc_qi.Rd: \value
  calc_sigma.Rd: \value
  calc_XHiY.Rd: \value
  center_kinship.Rd: \value
  eigen_proc.Rd: \value
  eigen2.Rd: \value
  MphCalcLogL.Rd: \value
  update_e.Rd: \value
  update_u.Rd: \value
  update_v.Rd: \value
  UpdateRL_B.Rd: \value

> Please add and explain the function's results.

> Please fix and resubmit.

> Best,
> Swetlana Herbrandt

