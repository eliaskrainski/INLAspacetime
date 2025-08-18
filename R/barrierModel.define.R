#' Define a spacetime model object for the `f()` call.
#'
#' @param mesh a spatial mesh
#' @param barrier.triangles a integer vector to specify which
#' triangles centers are in the barrier domain,
#' or a list with integer vector if more than one.
#' @param prior.range numeric vector containing U and a
#' to define the probability statements P(range < U) = a
#' used to setup the PC-prior for range.
#' If a = 0 or a = NA, then U is taken to be the fixed value for the range.
#' @param prior.sigma numeric vector containing U and a
#' to define the probability statements P(range > U) = a
#' used to setup the PC-prior for sigma.
#' If a = 0 or a = NA, then U is taken to be the fixed value for sigma.
#' @param range.fraction numeric to specify the fraction of the range
#' for the barrier domain. Default value is 0.1.
#' This has to be specified with care in order to have it small enough
#' to make it act as barrier but not too small in order to
#' prevent numerical issues.
#' @param constr logical, default is FALSE, to indicate if the
#' integral of the field over the domain is to be constrained to zero.
#' @param debug integer, default is zero, indicating the verbose level.
#' Will be used as logical by INLA.
#' @param useINLAprecomp logical, default is TRUE, indicating if it is to
#' be used the shared object pre-compiled by INLA.
#' This is not considered if 'libpath' is provided.
#' @param libpath string, default is NULL, with the path to the shared object.
#' @details
#' See the paper.
#' @return objects to be used in the f() formula term in INLA.
#' @export
barrierModel.define <-
  function(mesh, barrier.triangles,
           prior.range, prior.sigma,
           range.fraction = 0.1,
           constr = FALSE,
           debug = FALSE,
           useINLAprecomp = TRUE,
           libpath = NULL) {

    stopifnot(length(prior.range)==2)
    prior.range <- as.numeric(prior.range)
    if(is.na(prior.range[2])) {
      prior.range[2] <- 0
    }
    stopifnot(prior.range[1] > 0)
    stopifnot(prior.range[2] >= 0)
    stopifnot(prior.range[2] < 1)

    stopifnot(length(prior.sigma)==2)
    prior.sigma <- as.numeric(prior.sigma)
    if(is.na(prior.sigma[2])) {
      prior.sigma[2] <- 0
    }
    stopifnot(prior.sigma[1] > 0)
    stopifnot(prior.sigma[2] >= 0)
    stopifnot(prior.sigma[2] < 1)

    INLAversion <- INLAtools::packageCheck(
      name = "INLA",
      minimum_version = "24.10.07",
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
      hasverbose <- FALSE
    }
    stopifnot(file.exists(libpath))

    bfem <- INLAspacetime::mesh2fem.barrier(mesh, barrier.triangles)
    n <- nrow(bfem$I)

    if(!is.list(barrier.triangles)) {
      barrier.triangles <- list(barrier.triangles)
    }
    no <- length(barrier.triangles) + 1
    if(length(range.fraction) == 1) {
      range.fraction <- rep(range.fraction, no-1)
    } else {
      stopifnot(length(range.fraction)==(no-1))
    }

    Imat <- bfem$I
    Dmat <- bfem$D[[1]]
    CC <- bfem$C[[1]]
    for(o in 2:no) {
      CC <- CC + bfem$C[[o]] * (range.fraction[o-1]^2)
      Dmat <- Dmat +  bfem$D[[o]] * (range.fraction[o-1]^2)
    }
    iC <- Diagonal(n, 1 / CC)

    lmats <- INLAtools::upperPadding(
      list(
        ici = t(Imat) %*% iC %*% Imat,
        icd = t(Imat) %*% iC %*% Dmat,
        dci = t(Dmat) %*% iC %*% Imat,
        dcd = t(Dmat) %*% iC %*% Dmat
      ),
      relative = FALSE
    )
    stopifnot(n == nrow(lmats$graph))

    args0 <- list(
      model = "inla_cgeneric_barrier",
      shlib = libpath,
      n = as.integer(n),
      debug = as.integer(debug)
    )
    if(hasverbose) { ## to work with old C versions
      args0$verbose <- as.integer(0)
    }
    if(INLAversion<="25.02.10") {
      args0$prs <- prior.range
    } else {
      args0$prange <- prior.range
    }
    args0$psigma <- prior.sigma

    the_model <- do.call(
      "inla.cgeneric.define",
      c(args0,
        list(
        ii = lmats$graph@i,
        jj = lmats$graph@j,
        xx = t(lmats$xx)
      )
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


    the_model$mapper <- NULL
    if(requireNamespace("inlabru")) {
      if(!is.na(INLAtools::packageCheck(
        name = "inlabru",
        minimum_version = "2.12.0.9021",
        quietly = TRUE
      ))) {
        the_model$mapper <- inlabru::bru_mapper(mesh)
      }
    }

    return(the_model)
  }
