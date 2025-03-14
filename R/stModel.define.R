#' Define a spacetime model object for the `f()` call.
#'
#' @param smesh a spatial mesh
#' @param tmesh a temporal mesh
#' @param model a three characters string to specify the
#' smoothness alpha (each one as integer) parameters.
#' Currently it considers the `102`, `121`, `202` and `220` models.
#' @param control.priors a named list with parameter priors,
#' named as `prs`, `prt` and `psigma`, each one as a vector
#' with length two containing (U, a) to define the
#' corresponding PC-prior such that, respectively,
#' P(range.spatial<U)=a, P(range.temporal<U)=a or P(sigma>U)=a.
#' If a=0 then U is taken to be the fixed value of the parameter.
#' @param constr logical to indicate if the integral of the field
#' over the domain is to be constrained to zero. Default value is FALSE.
#' @param debug integer indicating the verbose level.
#' Will be used as logical by INLA.
#' @param useINLAprecomp logical indicating if is to be used
#' shared object pre-compiled by INLA. Not considered if
#' libpath is provided.
#' @param libpath string to the shared object. Default is NULL.
#' @details
#' The matrices of the kronecker product in Theorem 4.1 of
#' Lindgren et. al. (2024) are computed with the
#' [stModel.matrices] and the parameters are as Eq (19-21),
#' but parametrized in log scale.
#' @return objects to be used in the f() formula term in INLA.
#' @references
#' Finn Lindgren, Haakon Bakka, David Bolin, Elias Krainski and Håvard Rue (2024).
#' A diffusion-based spatio-temporal extension of Gaussian Matérn fields.
#' [SORT 48 (1), 3-66](https://www.idescat.cat/sort/sort481/48.1.1.Lindgren-etal.pdf)
#' <doi: 10.57645/20.8080.02.13>
#' @export
#' @importFrom fmesher fm_manifold
stModel.define <-
  function(smesh, tmesh, model, control.priors,
           constr = FALSE, debug = FALSE,
           useINLAprecomp = TRUE, libpath = NULL) {
    stopifnot(model %in% c("102", "121", "202", "220"))

    stopifnot(fm_manifold(smesh, c("S", "R")))
    Rmanifold <- fm_manifold(smesh, "R") + 0L

    dimension <- fm_manifold_dim(smesh)
    stopifnot(dimension > 0)

    stopifnot(length(control.priors$prt)==2)
    stopifnot(length(control.priors$prs)==2)
    stopifnot(length(control.priors$psigma)==2)

    stopifnot(control.priors$prt[1]>0)
    stopifnot(control.priors$prs[1]>0)
    stopifnot(control.priors$psigma[1]>0)

    stopifnot(control.priors$prt[2]>=0)
    stopifnot(control.priors$prs[2]>=0)
    stopifnot(control.priors$psigma[2]>=0)

    stopifnot(control.priors$prt[2]<1)
    stopifnot(control.priors$prs[2]<1)
    stopifnot(control.priors$psigma[2]<1)

    INLAversion <- check_package_version_and_load(
      pkg = "INLA",
      minimum_version = "23.08.16",
      quietly = TRUE
    )
    if (is.null(libpath)) {
      if(length(useINLAprecomp)>1) {
        warning("length(useINLAprecomp)>1, first taken!")
        useINLAprecomp <- useINLAprecomp[1]
      }
      stopifnot(is.logical(useINLAprecomp))
      if(is.na(INLAversion) & useINLAprecomp) {
        stop("Update INLA or try `useINLAprecomp = FALSE`!")
      }
      if (useINLAprecomp) {
        hasverbose <- (INLAversion<="25.02.10") ## to work with old C versions
        libpath <- INLA::inla.external.lib("INLAspacetime")
      } else {
        hasverbose <- FALSE
        libpath <- system.file("libs", package = "INLAspacetime")
        if (Sys.info()["sysname"] == "Windows") {
          libpath <- file.path(libpath, "INLAspacetime.dll")
        } else {
          libpath <- file.path(libpath, "INLAspacetime.so")
        }
      }
    } else {
      hasverbose <- FALSE ## assumed...
    }
    stopifnot(file.exists(libpath))

    alphas <- as.integer(strsplit(model, "")[[1]])
    alpha <- alphas[3] + alphas[2] * (alphas[1] - 0.5)
    nu.s <- alpha - dimension / 2
    nu.t <- min(alphas[1] - 0.5, nu.s / alphas[2])

    if (debug) {
      print(c(alphas = alphas))
      print(c(alpha = alpha, nu.s = nu.s, nu.t = nu.t))
    }

    log.C.t <- lgamma(alphas[1] - 0.5) - lgamma(alphas[1]) - 0.5 * log(4 * pi)
    cc <- c(
      c1 = 0.5 * log(8 * nu.s),
      c2 = -0.5 * log(8 * nu.t),
      c3 = NA
    )
    if (Rmanifold) {
      log.C.Rd <- lgamma(alpha - (dimension * 0.5)) - lgamma(alpha) -
        (dimension * 0.5) * log(4 * pi)
      cc[3] <- 0.5 * (log.C.t + log.C.Rd)
      if (debug) {
        cat("R manifold, cc[3] = ", cc[3], "\n")
      }
    } else {
      log.C.S2.part <- -log(4 * pi) ## S1???
      cc[3] <- 0.5 * (log.C.t + log.C.S2.part
      ) ## c3 part for S2 (S1???), to be completed in C
      if (debug) {
        cat("S manifold, cc[3] = ", cc[3], "\n")
      }
    }
    if (debug) {
      print(c(cc = cc))
    }

    mm <- stModel.matrices(smesh, tmesh, model, constr)
    n <- smesh$n * tmesh$n
    nm <- ncol(mm$TT)
    stopifnot(nm == length(mm$bb))
    jmm <- pmatch(paste0("M", 1:nm), names(mm))
    stopifnot(length(jmm[complete.cases(jmm)]) == nm)

    lmats <- upperPadding(mm[jmm], relative = FALSE)
    stopifnot(n == nrow(lmats$graph))

    args0 <- list(
      model = "inla_cgeneric_sstspde",
      shlib = libpath,
      n = as.integer(n),
      debug = as.integer(debug)
      )
    if(hasverbose) { ## to work with old C versions
      args0$verbose <- as.integer(0)
    }

    the_model <- do.call(
      "inla.cgeneric.define",
      c(args0,
        list(
        Rmanifold = as.integer(Rmanifold),
        dimension = as.integer(dimension),
        aaa = as.integer(alphas),
        nm = as.integer(nm),
        cc = as.double(cc),
        bb = mm$bb,
        prs = control.priors$prs,
        prt = control.priors$prt,
        psigma = control.priors$psigma,
        ii = lmats$graph@i,
        jj = lmats$graph@j,
        tt = t(mm$TT),
        xx = t(lmats$xx)
      ))
    )
    if (constr) {
      the_model$f$extraconstr <- mm$extraconstr
    }
    # Prepend specialised model class identifier, for bru_mapper use:
    class(the_model) <- c("stModel_cgeneric", class(the_model))
    # Add objects needed by bru_get_mapper.stModel_cgeneric:
    # (alternatively, construct the mapper already here, but that would
    # require loading inlabru even when it's not going to be used)
    the_model[["smesh"]] <- smesh
    the_model[["tmesh"]] <- tmesh

    the_model
  }
