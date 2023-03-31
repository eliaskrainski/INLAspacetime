
<!-- README.md is generated from README.Rmd. Please edit that file -->

# INLAspacetime

<!-- badges: start -->

[![CRAN
Status](http://www.r-pkg.org/badges/version-last-release/INLAspacetime)](https://cran.r-project.org/package=INLAspacetime)
[![INLAspacetime status
badge](https://eliaskrainski.r-universe.dev/badges/INLAspacetime)](https://eliaskrainski.r-universe.dev)
[![R build
status](https://github.com/eliaskrainski/INLAspacetime/workflows/R-CMD-check/badge.svg)](https://github.com/eliaskrainski/INLAspacetime/actions)
[![R code coverage
status](https://github.com/eliaskrainski/INLAspacetime/workflows/test-coverage/badge.svg)](https://github.com/eliaskrainski/INLAspacetime/actions)
[![lintr
status](https://github.com/eliaskrainski/INLAspacetime/workflows/lint/badge.svg)](https://github.com/eliaskrainski/INLAspacetime/actions)
[![Codecov test
coverage](https://codecov.io/gh/eliaskrainski/INLAspacetime/branch/devel/graph/badge.svg)](https://app.codecov.io/gh/eliaskrainski/INLAspacetime?branch=devel)
<!-- badges: end -->

This is a R package to implement certain spatial and spatio-temporal
models taking use to the `cgeneric` interface in the INLA package. This
interface is a way to implement models by writing `C` code to build the
precision matrix compiling it so that INLA can use it internally.

## We have implemented

1.  some of the models presented in <https://arxiv.org/abs/2006.04917>

2.  the barrier model proposed in
    <https://doi.org/10.1016/j.spasta.2019.01.002>

## Installation

<!-- You can install the current [CRAN](https://CRAN.R-project.org) version of INLAspacetime: -->
<!-- ```{r cran-installation, eval = FALSE} -->
<!-- install.packages("INLAspacetime") -->
<!-- ``` -->

You can install the latest version of INLAspacetime from
[GitHub](https://github.com/eliaskrainski/INLAspacetime) with

``` r
## install.packages("remotes")
remotes::install_github("eliaskrainski/INLAspacetime")
```

<!-- or track the development version builds via [inlabru-org.r-universe.dev](https://inlabru-org.r-universe.dev/ui#builds): -->
<!-- ```{r universe-installation, eval = FALSE} -->
<!-- ## Enable universe(s) by inlabru-org -->
<!-- options(repos = c( -->
<!--   inlabruorg = "https://inlabru-org.r-universe.dev", -->
<!--   INLA = "https://inla.r-inla-download.org/R/testing", -->
<!--   CRAN = "https://cloud.r-project.org" -->
<!-- )) -->
<!-- ## Install it -->
<!-- install.packages("INLAspacetime") -->
<!-- ``` -->

# See the vignettes for examples

We will have tutorials and examples at
<https://eliaskrainski.github.io/INLAspacetime/>

# Example

This is a basic example which fit a spacetime model for some fake data.
The model fitting using **inlabru** facilitates coding.

``` r
n <- 5
dataf <- data.frame(
    s1   = runif(n, -1, 1),
    s2   = runif(n, -1, 1),
    time = runif(n, 1, 4),
    y    = rnorm(n, 0, 1))
str(dataf)
#> 'data.frame':    5 obs. of  4 variables:
#>  $ s1  : num  0.588 -0.453 -0.973 0.398 0.833
#>  $ s2  : num  -0.7447 -0.5511 -0.0729 0.4581 0.9561
#>  $ time: num  3.62 1.42 2.55 3.99 1.25
#>  $ y   : num  0.157 1.385 2.036 -0.295 0.963
```

Loading the packages:

``` r
library(INLA)
#> Loading required package: Matrix
#> Loading required package: foreach
#> Loading required package: parallel
#> Loading required package: sp
#> This is INLA_23.03.26 built 2023-03-26 19:17:37 UTC.
#>  - See www.r-inla.org/contact-us for how to get help.
#>  - To enable PARDISO sparse library; see inla.pardiso()
library(INLAspacetime)
library(inlabru)
```

Define spatial and temporal discretization meshes

``` r
smesh <- inla.mesh.2d(
  loc = cbind(0,0), 
  max.edge = 5, 
  offset = 2)
tmesh <- inla.mesh.1d(
  loc = 0:5)
```

Define the spacetime model object to be used

``` r
stmodel <- stModel.define(
    smesh = smesh, ## spatial mesh
    tmesh = tmesh, ## temporal mesh
    model = '121', ## model, see the paper
    control.priors = list(
        prs = c(1, 0.1), ## P(spatial range < 1) = 0.1
        prt = c(5, 0), ## fixed to 5
        psigma = c(1, 0.1) ## P(sigma > 1) = 0.1
        )
    )
```

Define the data model: the linear predictor terms

``` r
linpred <- ~ 1 +
    field(list(space = cbind(s1, s2), 
               time = time),
          model = stmodel)
```

Setting the likelihood

``` r
likeprec <- list(prec = list(
  initial = 10, fixed = FALSE))
datalike <- like(
  formula = y ~ ., 
  family = "gaussian",
  control.family = list(
    hyper = likeprec), ## TO DO: not going through bru
  data=dataf)
```

Fitting

``` r
result <- 
  bru(
    components = linpred,
    datalike,
    options = list(
      control.inla = list(
        int.strategy = "eb"
        ),
      verbose = !TRUE)
    )
#> Warning in inla.model.properties.generic(inla.trim.family(model), mm[names(mm) == : Model 'cgeneric' in section 'latent' is marked as 'experimental'; changes may appear at any time.
#>   Use this model with extra care!!! Further warnings are disabled.
```

Summary of the model parameters

``` r
result$summary.fixed
#>                mean       sd 0.025quant  0.5quant 0.975quant      mode kld
#> Intercept 0.5619682 2.149405  -3.650789 0.5619682   4.774725 0.5619682   0
result$summary.hyperpar
#>                                                  mean           sd   0.025quant
#> Precision for the Gaussian observations 18329.3225260 1.807974e+04 1233.6322877
#> Theta1 for field                            2.3183849 3.508814e-01    1.5331198
#> Theta2 for field                           -0.3303781 3.479795e-01   -0.8930268
#>                                              0.5quant   0.975quant        mode
#> Precision for the Gaussian observations 12865.9846347 6.651524e+04 3388.114014
#> Theta1 for field                            2.3553357 2.890631e+00    2.515316
#> Theta2 for field                           -0.3689591 4.502314e-01   -0.534043
```
