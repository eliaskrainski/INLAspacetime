
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
models, including some of the spatio-temporal models proposed
[here](https://www.idescat.cat/sort/sort481/48.1.1.Lindgren-etal.pdf).
It uses the `cgeneric` interface in the INLA package, to implement
models by writing `C` code to build the precision matrix compiling it so
that INLA can use it internally.

## We have implemented

1.  some of the models presented in *A diffusion-based spatio-temporal
    extension of Gaussian Matérn fields* (2024). *Finn Lindgren, Haakon
    Bakka, David Bolin, Elias Krainski and Håvard Rue*. SORT 48 (1)
    January-June 2024, 3-66.
    (<https://www.idescat.cat/sort/sort481/48.1.1.Lindgren-etal.pdf>)

2.  the barrier (and transparent barriers) model proposed in
    <https://doi.org/10.1016/j.spasta.2019.01.002>

## Vignettes

Please check
[here](https://eliaskrainski.github.io/INLAspacetime/articles/)

## Installation

The ‘INLA’ package is a suggested one, but you will need it for actually
fitting a model. You can install it with

``` r
install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/testing"), dep=TRUE) 
```

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

# A spacetime example

Simulate some fake data.

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

Loading packages:

``` r
library(fmesher)
library(INLA)
library(INLAspacetime)
#> see more on https://eliaskrainski.github.io/INLAspacetime
```

Define spatial and temporal discretization meshes

``` r
smesh <- fm_mesh_2d(
  loc = cbind(0,0), 
  max.edge = 5, 
  offset = 2)
tmesh <- fm_mesh_1d(
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
#> Warning in stModel.define(smesh = smesh, tmesh = tmesh, model = "121",
#> control.priors = list(prs = c(1, : Setting 'useINLAprecomp = FALSE' to use new
#> code.
```

## Fit the model

Define a projector matrix from the spatial and temporal meshes to the
data

``` r
Aproj <- inla.spde.make.A(
    mesh = smesh,
    loc = cbind(dataf$s1, dataf$s2),
    group = dataf$time,
    group.mesh = tmesh
)
```

Create a ‘fake’ column to be used as index. in the `f()` term

``` r
dataf$st <- NA
```

Setting the likelihood precision (as fixed)

``` r
ctrl.lik <- list(
  hyper = list(
    prec = list(
      initial = 10, 
      fixed = TRUE)    
  )
)
```

Combine a ‘fake’ index column with `A.local`

``` r
fmodel <- y ~ f(st, model = stmodel, A.local = Aproj)
```

Call the main `INLA` function:

``` r
fit <- inla(
    formula = fmodel,
    data = dataf,
    control.family = ctrl.lik)
```

Posterior marginal summaries for fixed effect and the model parameters
that were not fixed.

``` r
fit$summary.fixed
#>                 mean      sd 0.025quant  0.5quant 0.975quant      mode
#> (Intercept) 0.693389 4.03265  -6.962331 0.5227188   9.417425 0.5550712
#>                      kld
#> (Intercept) 7.398472e-05
fit$summary.hyperpar
#>                   mean        sd 0.025quant 0.5quant 0.975quant     mode
#> Theta1 for st 1.199222 0.4918533  0.3653818 1.161539   2.277396 0.974993
#> Theta2 for st 1.435517 0.1710676  1.1031120 1.434032   1.776667 1.427752
```

## Using the **inlabru**

``` r
library(inlabru)
```

Setting the observation (likelihood) model object

``` r
data_model <- bru_obs(
  formula = y ~ ., 
  family = "gaussian",
  control.family = ctrl.lik, 
  data = dataf)
```

Define the data model: the linear predictor terms

``` r
linpred <- ~ 1 +
    field(list(space = cbind(s1, s2), 
               time = time),
          model = stmodel)
```

Fitting

``` r
result <- bru(
  components = linpred,
  data_model)
```

Summary of the model parameters

``` r
result$summary.fixed
#>                mean       sd 0.025quant 0.5quant 0.975quant      mode
#> Intercept 0.6690302 3.970182  -6.887199 0.509471   9.214066 0.5379221
#>                    kld
#> Intercept 5.683968e-05
result$summary.hyperpar
#>                      mean        sd 0.025quant 0.5quant 0.975quant      mode
#> Theta1 for field 1.190438 0.4868951  0.3623876 1.153809   2.256071 0.9726162
#> Theta2 for field 1.435268 0.1709839  1.1033563 1.433674   1.776580 1.4269195
```
