#' Define a spacetime model object for the `f()` call.
#'
#' @param mesh a spatial mesh
#' @param barrier.triangles a integer vector to specify which
#' triangles centers are in the barrier domain
#' @param prior.range numeric vector containing U and a
#' to define the probability statements P(range < U) = a
#' used to setup the PC-prior for range.
#' If a = 0 then U is taken to be the fixed value for the range.
#' @param prior.sigma numeric vector containing U and a
#' to define the probability statements P(range > U) = a
#' used to setup the PC-prior for sigma.
#' If a = 0 then U is taken to be the fixed value for sigma.
#' @param range.fraction numeric to specify the fraction of the range
#' for the barrier domain. Default value is 0.2.
#' @param constr logical to indicate if the integral of the field
#' over the domain is to be constrained to zero. Default value is FALSE.
#' @param debug logical indicating if to run in debug mode.
#' @param verbose logical indicating if to print parameter values.
#' @param useINLAprecomp logical indicating if is to be used
#' shared object pre-compiled by INLA. Not considered if
#' libpath is provided.
#' @param libpath string to the shared object. Default is NULL.
#' @details
#' See the paper.
#' @return objects to be used in the f() formula term in INLA.
#' @export
barrierModel.define <-
  function(mesh, barrier.triangles,
           prior.range, prior.sigma, range.fraction = 0.2,
           constr = FALSE, debug = FALSE, verbose = FALSE,
           useINLAprecomp = TRUE, libpath = NULL) {
    stopifnot(all(c(
      length(prior.range),
      length(prior.sigma) > 1
    )))
    stopifnot(prior.range[1] > 0)
    stopifnot(prior.range[2] >= 0)
    stopifnot(prior.range[2] < 1)
    stopifnot(prior.sigma[2] >= 0)
    stopifnot(prior.sigma[2] < 1)
    stopifnot(prior.sigma[1] > 0)

    bfem <- INLAspacetime::mesh2fem.barrier(mesh, barrier.triangles)
    n <- nrow(bfem$I)

    Imat <- bfem$I
    Dmat <- bfem$D[[1]] + bfem$D[[2]] * range.fraction^2
    iC <- Diagonal(n, 1 / (bfem$C[[1]] + bfem$C[[2]] * range.fraction^2))

    lmats <- upperPadding(
      list(
        ici = t(Imat) %*% iC %*% Imat,
        icd = t(Imat) %*% iC %*% Dmat,
        dci = t(Dmat) %*% iC %*% Imat,
        dcd = t(Dmat) %*% iC %*% Dmat
      ),
      relative = FALSE
    )
    stopifnot(n == nrow(lmats$graph))

    if (is.null(libpath)) {
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

    the_model <- do.call(
      "inla.cgeneric.define",
      list(
        model = "inla_cgeneric_barrier",
        shlib = libpath,
        n = n,
        debug = as.integer(debug),
        verbose = as.integer(verbose),
        prs = prior.range, ## TO DO: name on current CRAN version, change to 'prange'
        psigma = prior.sigma,
        ii = lmats$graph@i,
        jj = lmats$graph@j,
        xx = t(lmats$xx)
      )
    )
    if (constr) {
      the_model$f$extraconstr <- list(
        A = matrix(1 / n, 1, n), e = 0.0
      )
    }
    # Prepend specialised model class identifier, for bru_mapper use:
    class(the_model) <- c("barrierModel_cgeneric", class(the_model))
    # Add objects needed by bru_get_mapper.barrierModel_cgeneric:
    # (alternatively, construct the mapper already here, but that would
    # require loading inlabru even when it's not going to be used)
    the_model[["mesh"]] <- mesh

    the_model
  }
