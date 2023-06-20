
<!-- README.md is generated from README.Rmd. Please edit that file -->

# INLAspacetime

<!-- badges: start -->

[![CRAN
Status](http://www.r-pkg.org/badges/version-last-release/INLAspacetime)](https://cran.r-project.org/package=INLAspacetime)
[![](https://cranlogs.r-pkg.org/badges/INLAspacetime)](https://cran.r-project.org/package=INLAspacetime)
[![check
no-suggestions](https://github.com/eliaskrainski/INLAspacetime/workflows/R-CMD-check-no-suggests/badge.svg)](https://github.com/eliaskrainski/INLAspacetime/actions)
[![check](https://github.com/eliaskrainski/INLAspacetime/workflows/R-CMD-check/badge.svg)](https://github.com/eliaskrainski/INLAspacetime/actions)
[![pkgdown](https://github.com/eliaskrainski/INLAspacetime/workflows/pkgdown/badge.svg)](https://github.com/eliaskrainski/INLAspacetime/actions)
<!-- badges: end -->

This is a R package to implement certain spatial and spatio-temporal
models taking use to the `cgeneric` interface in the INLA package. This
interface is a way to implement models by writing `C` code to build the
precision matrix compiling it so that INLA can use it internally.

## Installation

You can install the current [CRAN](https://CRAN.R-project.org) version
of INLAspacetime:

``` r
install.packages("INLAspacetime")
```

You can install the latest version of INLAspacetime from
[GitHub](https://github.com/eliaskrainski/INLAspacetime) with

``` r
## install.packages("remotes")
remotes::install_github("eliaskrainski/INLAspacetime",  build_vignettes=TRUE)
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

## We have implemented

1.  some of the models presented in <https://arxiv.org/abs/2006.04917>

2.  the barrier model proposed in
    <https://doi.org/10.1016/j.spasta.2019.01.002>

# Example

This is a basic example which fit a spacetime model for some fake data.
The model fitting using **inlabru** facilitates coding.

``` r
set.seed(1)
n <- 5
dataf <- data.frame(
    s1   = runif(n, -1, 1),
    s2   = runif(n, -1, 1),
    time = runif(n, 1, 4),
    y    = rnorm(n, 0, 1))
str(dataf)
#> 'data.frame':    5 obs. of  4 variables:
#>  $ s1  : num  -0.469 -0.256 0.146 0.816 -0.597
#>  $ s2  : num  0.797 0.889 0.322 0.258 -0.876
#>  $ time: num  1.62 1.53 3.06 2.15 3.31
#>  $ y   : num  -0.00577 2.40465 0.76359 -0.79901 -1.14766
```

Loading the packages:

``` r
library(INLA)
#> Loading required package: Matrix
#> Loading required package: foreach
#> Loading required package: parallel
#> Loading required package: sp
#> The legacy packages maptools, rgdal, and rgeos, underpinning this package
#> will retire shortly. Please refer to R-spatial evolution reports on
#> https://r-spatial.org/r/2023/05/15/evolution4.html for details.
#> This package is now running under evolution status 0
#> This is INLA_23.05.30-1 built 2023-05-30 11:52:19 UTC.
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
        prt = c(5, 0), ## temporal range fixed to 5
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
ctrlf <- list(
  hyper = list(
    prec = list(
      initial = 10, 
      fixed = TRUE)    
  )
)
datalike <- like(
  formula = y ~ ., 
  family = "gaussian",
  control.family = ctrlf, 
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
#>                mean      sd 0.025quant  0.5quant 0.975quant      mode kld
#> Intercept 0.5211454 3.49482  -6.328575 0.5211454   7.370866 0.5211454   0
result$summary.hyperpar
#>                      mean        sd 0.025quant 0.5quant 0.975quant      mode
#> Theta1 for field 1.202935 0.5125858  0.3321324 1.162519   2.324768 0.9707531
#> Theta2 for field 1.435746 0.1706775  1.1044643 1.434143   1.776487 1.4273398
```

## Vignettes

Please check it out at the
[Tutorials](https://eliaskrainski.github.io/INLAspacetime/)
