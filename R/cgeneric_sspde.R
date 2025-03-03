#' Define the stationary SPDE cgeneric model for INLA.
#' @param mesh triangulation mesh to discretize the model.
#' @param alpha integer used to compute the smoothness parameter.
#' @param control.priors named list with parameter priors.
#' This shall contain prior.range and prior.sigma as length two
#' vectors with (U, a) to define the PC-prior parameters such that
#' P(range<U)=a and P(sigma>U)=a, respectively.
#' See Fuglstad et. al. (2019) <DOI: 10.1080/01621459.2017.1415907>.
#' If a=0 then U is taken to be the fixed value of the parameter.
#' @param constr logical to indicate if the integral of the field
#' over the domain is to be constrained to zero. Default value is FALSE.
#' @param debug integer indicating the debug level.
#' Will be used as logical by INLA.
#' @param useINLAprecomp logical indicating if is to be used
#' shared object pre-compiled by INLA. Not considered if
#' libpath is provided.
#' @param libpath string to the shared object. Default is NULL.
#' @note
#' This is the stationary case of [INLA::inla.spde2.pcmatern()]
#' with slight change on the marginal variance when the domain is
#' the sphere, following Eq. (23) in Lindgren et. al. (2024).
#' @references
#' Geir-Arne Fuglstad, Daniel Simpson, Finn Lindgren & Håvard Rue (2019).
#' Constructing Priors that Penalize the Complexity of Gaussian Random Fields.
#' Journal of the American Statistical Association, V. 114, Issue 525.
#'
#' Finn Lindgren, Haakon Bakka, David Bolin, Elias Krainski and Håvard Rue (2024).
#' A diffusion-based spatio-temporal extension of Gaussian Matérn fields.
#' [SORT 48 (1), 3-66](https://www.idescat.cat/sort/sort481/48.1.1.Lindgren-etal.pdf)
#' <doi: 10.57645/20.8080.02.13>
#' @return objects to be used in the f() formula term in INLA.
#' @importFrom fmesher fm_manifold
#' @export
cgeneric_sspde <-
  function(mesh,
           alpha,
           control.priors,
           constr = FALSE,
           debug = FALSE,
           useINLAprecomp = TRUE,
           libpath = NULL) {

    stopifnot(fm_manifold(mesh, c("S", "R")))
    Rmanifold <- fm_manifold(mesh, "R") + 0L

    dimension <- fm_manifold_dim(mesh)
    stopifnot(dimension > 0)

    if (is.null(libpath)) {
      pkgs <- installed.packages()
      iINLAi <- which(pkgs[, "Package"] == "INLA")
      if(length(iINLAi)==0)
        stop("You need to install `INLA`!")
      INLAVersion <- pkgs[iINLAi, "Version"]
      if(useINLAprecomp & (!(INLAVersion>"25.02.10"))) {
        warning("Setting 'useINLAprecomp = FALSE' to use new code.")
        useINLAprecomp <- FALSE
      }
      if (useINLAprecomp) {
        libpath <- INLA::inla.external.lib("INLAspacetime")
      } else {
        libpath <- system.file("libs", package = "INLAspacetime")
        if (Sys.info()["sysname"] == "Windows") {
          libpath <- file.path(libpath, "INLAspacetime.dll")
        } else {
          libpath <- file.path(libpath, "INLAspacetime.so")
        }
      }
    }
    stopifnot(file.exists(libpath))

    if(alpha!=2) stop("Only 'alpha = 2' is supported for now.")
    nu.s <- alpha - dimension / 2

    if (debug) {
      print(c(alpha = alpha, nu.s = nu.s))
    }

    ## user parameters are {range, sigma}
    ## internal parameters are log(range, sigma)
    ## SPDE parameters {kappa, tau}:
    ##  \tau(\kappa^2 - D2)^{\alpha/2}u(s) = W(s)

    cc <- c(
      c1 = log(8 * nu.s),
      c2 = NA
    )
    if (Rmanifold) {
      cc[2] <- lgamma(alpha - (dimension * 0.5)) -
        lgamma(alpha) -(dimension * 0.5) * log(4 * pi)
      if (debug) {
        cat("R manifold, cc[2] = ", cc[2], "\n")
      }
    } else {
      ## c2 part depending on gamma to be completed in C
      cc[2] <- -log(4 * pi) ## (S1???)
      if (debug) {
        cat("S manifold, cc[2] = ", cc[2], "\n")
      }
    }
    if (debug) {
      print(c(cc = cc))
    }

    fem <- fmesher::fm_fem(mesh, order = 2)
    stopifnot((n <- nrow(fem$g1))>0)
    lmats <- upperPadding(
      fem[c("c0", "g1", "g2")],
      relative = FALSE)
    stopifnot(n == nrow(lmats$graph))

    the_model <- do.call(
      "inla.cgeneric.define",
      list(
        model = "inla_cgeneric_sspde",
        shlib = libpath,
        n = as.integer(n),
        debug = as.integer(debug),
        Rmanifold = as.integer(Rmanifold),
        dimension = as.integer(dimension),
        alpha = as.integer(alpha),
        nm = as.integer(ncol(lmats$xx)),
        ii = as.integer(lmats$graph@i),
        jj = as.integer(lmats$graph@j),
        cc = as.double(cc),
        prange = as.double(control.priors$prange),
        psigma = as.double(control.priors$psigma),
        xx = t(lmats$xx)
        )
    )
    if (constr) {
      the_model$f$extraconstr <- list(
        A = matrix(fem$c0@x, nrow = 1), e = 0)
    }
    the_model$"mesh" <- mesh
    the_model$"fem" <- fem

    return(the_model)
  }
