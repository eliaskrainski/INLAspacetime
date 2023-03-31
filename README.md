
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
<https://eliaskrainski.github.io/INLAspacetime/articles/>

# Example

This is a basic example which fit a spacetime model for some fake data.
The model fitting using **inlabru** facilitates coding.

``` r
n <- 5
dataf <- data.frame(
    s1=runif(n, -1, 1),
    s2=runif(n, -1, 1),
    time=runif(n, 1, 4),
    y=rnorm(n, 0, 1))
str(dataf)
#> 'data.frame':    5 obs. of  4 variables:
#>  $ s1  : num  -0.571 0.447 -0.442 -0.386 -0.63
#>  $ s2  : num  -0.937 -0.194 0.267 -0.285 0.947
#>  $ time: num  1.2 3.53 3.18 3.34 3.52
#>  $ y   : num  0.136 -0.17 -1.593 -0.794 -2.121
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
smesh <- inla.mesh.2d(cbind(0,0), max.edge=5, offset=2)
tmesh <- inla.mesh.1d(0:5)
```

Define the spacetime model object to be used

``` r
stmodel <- stModel.define(
    smesh, tmesh, '121', 
    control.priors=list(
        prs=c(1, 0.5),
        prt=c(5, 0.5),
        psigma=c(1, 0.5)))
```

Define the data model: the linear predictor terms

``` r
M <- ~ -1 + Intercept(1) +
    field(list(space = cbind(s1, s2), time=time),
          model=stmodel)
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
    components = M,
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
#>                 mean        sd 0.025quant   0.5quant 0.975quant       mode kld
#> Intercept -0.6008767 0.7269845   -2.02574 -0.6008767  0.8239867 -0.6008767   0
result$summary.hyperpar
#>                                                  mean           sd   0.025quant
#> Precision for the Gaussian observations  1.806235e+04 1.809052e+04 1232.7375106
#> Theta1 for field                         8.890266e-02 4.414888e-01   -0.8970600
#> Theta2 for field                        -3.945975e-01 1.100128e+00   -2.6341394
#> Theta3 for field                         1.154106e+00 3.986351e-01    0.4281037
#>                                              0.5quant   0.975quant         mode
#> Precision for the Gaussian observations 12582.6392583 6.637983e+04 3393.4055569
#> Theta1 for field                            0.1210754 8.377859e-01    0.2947973
#> Theta2 for field                           -0.3788068 1.741826e+00   -0.3108136
#> Theta3 for field                            1.1357565 2.005462e+00    1.0465864
```
