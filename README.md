
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
models taking use to the `cgeneric` interface in the INLA package,
including some of the spatio-temporal models proposed
[here](https://www.idescat.cat/sort/sort481/48.1.1.Lindgren-etal.pdf).
This interface is a way to implement models by writing `C` code to build
the precision matrix compiling it so that INLA can use it internally.

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

## We have implemented

1.  some of the models presented in
    <https://www.idescat.cat/sort/sort481/48.1.1.Lindgren-etal.pdf>

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
library(INLAspacetime)
#> Loading required package: fmesher
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
#> 
#>  *** inla.core.safe:  The inla program failed, but will rerun in case better initial values may help. try=1/1
#> Warning in iinla(model = info[["model"]], lhoods = info[["lhoods"]], options = info[["options"]]): iinla: Problem in inla: Error in inla.core.safe(formula = formula, family = family, contrasts = contrasts,  : 
#>   The inla-program exited with an error. Unless you interupted it yourself, please rerun with verbose=TRUE and check the output carefully.
#>   If this does not help, please contact the developers at <help@r-inla.org>.
#> The inla program failed and the maximum number of tries has been reached.
#> iinla: Problem in inla: 1: rmarkdown::render("/home/eliask/github/INLAspacetime/README.Rmd", 
#>        encoding = "UTF-8")
#> 2: knitr::knit(knit_input, knit_output, envir = envir, quiet = quiet)
#> 3: process_file(text, output)
#> 4: xfun:::handle_error(withCallingHandlers(if (tangle) process_tangle(group)  [...]
#>        error = function(e) if (xfun::pkg_available("rlang", "1.0.0")) rlang:: [...]
#>        function(loc) {
#>            setwd(wd)
#>            write_utf8(res, output %n% stdout())
#>            paste0("\nQuitting from lines ", loc)
#>        }, if (labels[i] != "") sprintf(" [%s]", labels[i]), get_loc)
#> 5: withCallingHandlers(if (tangle) process_tangle(group) else process_group(g [...]
#>        error = function(e) if (xfun::pkg_available("rlang", "1.0.0")) rlang:: [...]
#> 6: process_group(group)
#> 7: call_block(x)
#> 8: block_exec(params)
#> 9: eng_r(options)
#> 10: in_input_dir(evaluate(code, envir = env, new_device = FALSE, 
#>        keep_warning = if (is.numeric(options$warning)) TRUE else options$warning, 
#>        keep_message = if (is.numeric(options$message)) TRUE else options$message, 
#>        stop_on_error = if (is.numeric(options$error)) options$error else {
#>            if (options$error && options$include) 
#>                0L
#>            else 2L
#>        }, output_handler = knit_handlers(options$render, options)))
#> 11: in_dir(input_dir(), expr)
#> 12: evaluate(code, envir = env, new_device = FALSE, keep_warning = if (is.nume [...]
#>        keep_message = if (is.numeric(options$message)) TRUE else options$message, 
#>        stop_on_error = if (is.numeric(options$error)) options$error else {
#>            if (options$error && options$include) 
#>                0L
#>            else 2L
#>        }, output_handler = knit_handlers(options$render, options))
#> 13: evaluate::evaluate(...)
#> 14: withRestarts(with_handlers({
#>        for (expr in tle$exprs) {
#>            ev <- withVisible(eval(expr, envir))
#>            watcher$capture_plot_and_output()
#>            watcher$print_value(ev$value, ev$visible, envir)
#>        }
#>        TRUE
#>    }, handlers), eval_continue = function() TRUE, eval_stop = function() FALSE)
#> 15: withRestartList(expr, restarts)
#> 16: withOneRestart(withRestartList(expr, restarts[-nr]), restarts[[nr]])
#> 17: doWithOneRestart(return(expr), restart)
#> 18: withRestartList(expr, restarts[-nr])
#> 19: withOneRestart(expr, restarts[[1L]])
#> 20: doWithOneRestart(return(expr), restart)
#> 21: with_handlers({
#>        for (expr in tle$exprs) {
#>            ev <- withVisible(eval(expr, envir))
#>            watcher$capture_plot_and_output()
#>            watcher$print_value(ev$value, ev$visible, envir)
#>        }
#>        TRUE
#>    }, handlers)
#> 22: eval(call)
#> 23: eval(call)
#> 24: withCallingHandlers(code, message = function (cnd) 
#>    {
#>        watcher$capture_plot_and_output()
#>        if (on_message$capture) {
#>            watcher$push(cnd)
#>        }
#>        if (on_message$silence) {
#>            invokeRestart("muffleMessage")
#>        }
#>    }, warning = function (cnd) 
#>    {
#>        if (getOption("warn") >= 2 || getOption("warn") < 0) {
#>            return()
#>        }
#>        watcher$capture_plot_and_output()
#>        if (on_warning$capture) {
#>            cnd <- sanitize_call(cnd)
#>            watcher$push(cnd)
#>        }
#>        if (on_warning$silence) {
#>            invokeRestart("muffleWarning")
#>        }
#>    }, error = function (cnd) 
#>    {
#>        watcher$capture_plot_and_output()
#>        cnd <- sanitize_call(cnd)
#>        watcher$push(cnd)
#>        switch(on_error, continue = invokeRestart("eval_continue"), 
#>            stop = invokeRestart("eval_stop"), error = NULL)
#>    })
#> 25: withVisible(eval(expr, envir))
#> 26: eval(expr, envir)
#> 27: eval(expr, envir)
#> 28: bru(components = linpred, datalike, options = list(control.inla = list(int [...]
#>        verbose = !TRUE))
#> 29: iinla(model = info[["model"]], lhoods = info[["lhoods"]], options = info[[ [...]
#> 30: fm_try_callstack(...)
#> 31: do.call(INLA::inla, inla.options.merged, envir = environment(model$effects))
#> 32: (function (formula = NULL, family = "gaussian", contrasts = NULL, 
#>        data = NULL, quantiles = c(0.025, 0.5, 0.975), E = NULL, 
#>        offset = NULL, scale = NULL, weights = NULL, Ntrials = NULL, 
#>        strata = NULL, lp.scale = NULL, link.covariates = NULL, verbose = inla [...]
#>        lincomb = NULL, selection = NULL, control.compute = list(), 
#>        control.predictor = list(), control.family = list(), control.inla = list(), 
#>        control.fixed = list(), control.mode = list(), control.expert = list(), 
#>        control.hazard = list(), control.lincomb = list(), control.update = list(), 
#>        control.lp.scale = list(), control.pardiso = list(), only.hyperparam = [...]
#>        inla.call = inla.getOption("inla.call"), inla.arg = inla.getOption("in [...]
#>        num.threads = inla.getOption("num.threads"), keep = inla.getOption("keep"), 
#>        working.directory = inla.getOption("working.directory"), 
#>        silent = inla.getOption("silent"), inla.mode = inla.getOption("inla.mode"), 
#>        safe = inla.getOption("safe"), debug = inla.getOption("debug"), 
#>        .parent.frame = environment(formula)) 
#>    {
#>        set.warn <- function(a, b) {
#>            if (length(grep("(expand|E|weights|scale|offset|lp[.]scale)([.][.] [...]
#>                b)) == 0) {
#>                warning(paste0("Argument '", a, "=", b, "' expanded to NULL or [...]
#>                    "  This might be an error and you are requested to check t [...]
#>                    "  Move on with default values...\n"), immediate. = TRUE)
#>            }
#>            return(invisible())
#>        }
#>        is.set <- !is.null(substitute(E))
#>        nm <- if (is.set) 
#>            paste(collapse = " ", as.character(substitute(E)))
#>        else ""
#>        tmp <- try(eval(substitute(E), envir = data, enclos = .parent.frame), 
#>            silent = TRUE)
#>        E <- if (inherits(tmp, "try-error")) 
#>            NULL
#>        else tmp
#>        if (is.set && inherits(tmp, "try-error")) 
#>            set.warn("E", nm)
#>        is.set <- !is.null(substitute(offset))
#>        nm <- if (is.set) 
#>            paste(collapse = " ", as.character(substitute(offset)))
#>        else ""
#>        tmp <- try(eval(substitute(offset), envir = data, enclos = .parent.frame), 
#>            silent = TRUE)
#>        offset <- if (inherits(tmp, "try-error")) 
#>            NULL
#>        else tmp
#>        if (is.set && inherits(tmp, "try-error")) 
#>            set.warn("offset", nm)
#>        is.set <- !is.null(substitute(scale))
#>        nm <- if (is.set) 
#>            paste(collapse = " ", as.character(substitute(scale)))
#>        else ""
#>        tmp <- try(eval(substitute(scale), envir = data, enclos = .parent.frame), 
#>            silent = TRUE)
#>        scale <- if (inherits(tmp, "try-error")) 
#>            NULL
#>        else tmp
#>        if (is.set && inherits(tmp, "try-error")) 
#>            set.warn("scale", nm)
#>        is.set <- !is.null(substitute(weights))
#>        nm <- if (is.set) 
#>            paste(collapse = " ", as.character(substitute(weights)))
#>        else ""
#>        tmp <- try(eval(substitute(weights), envir = data, enclos = .parent.frame), 
#>            silent = TRUE)
#>        weights <- if (inherits(tmp, "try-error")) 
#>            NULL
#>        else tmp
#>        if (is.set && inherits(tmp, "try-error")) 
#>            set.warn("weights", nm)
#>        is.set <- !is.null(substitute(Ntrials))
#>        nm <- if (is.set) 
#>            paste(collapse = " ", as.character(substitute(Ntrials)))
#>        else ""
#>        tmp <- try(eval(substitute(Ntrials), envir = data, enclos = .parent.frame), 
#>            silent = TRUE)
#>        Ntrials <- if (inherits(tmp, "try-error")) 
#>            NULL
#>        else tmp
#>        if (is.set && inherits(tmp, "try-error")) 
#>            set.warn("Ntrials", nm)
#>        is.set <- !is.null(substitute(strata))
#>        nm <- if (is.set) 
#>            paste(collapse = " ", as.character(substitute(strata)))
#>        else ""
#>        tmp <- try(eval(substitut
#> iinla: Giving up and returning last successfully obtained result for diagnostic purposes.
```

Summary of the model parameters

``` r
result$summary.fixed
#> NULL
result$summary.hyperpar
#> NULL
```

## Vignettes

Please check it out at the
[Tutorials](https://eliaskrainski.github.io/INLAspacetime/)
