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
#' @param ... additional arguments, such as debug,useINLAprecomp,shlib.
#' @details
#' See the paper.
#' @return objects to be used in the f() formula term in INLA.
#' @export
barrierModel.define <-
  function(mesh, barrier.triangles,
           prior.range, prior.sigma,
           range.fraction = 0.1,
           constr = FALSE,
           ...) {

    dotArgs <- list(...)
    if(is.null(dotArgs$debug)) {
      dotArgs$debug <- FALSE
    } else {
      dotArgs$debug <- dotArgs$debug[1]
    }

    stopifnot(length(prior.range)==2)
    prior.range <- as.numeric(prior.range)
    if(is.na(prior.range[2])) {
      prior.range[2] <- 0
    }
    stopifnot(prior.range[1] > 0)
    stopifnot(prior.range[2] >= 0)
    stopifnot(prior.range[2] < 1)

    if(dotArgs$debug) {
      cat("The range prior parameters: ", prior.range, "\n")
    }

    stopifnot(length(prior.sigma)==2)
    prior.sigma <- as.numeric(prior.sigma)
    if(is.na(prior.sigma[2])) {
      prior.sigma[2] <- 0
    }
    stopifnot(prior.sigma[1] > 0)
    stopifnot(prior.sigma[2] >= 0)
    stopifnot(prior.sigma[2] < 1)

    if(dotArgs$debug) {
      cat("The sigma prior parameters: ", prior.sigma, "\n")
    }

    INLAversion <- packageCheck(
      name = "INLA",
      minimum_version = "24.10.07",
      quietly = TRUE
    )
    if(is.null(dotArgs$useINLAprecomp))
      dotArgs$useINLAprecomp <- TRUE
    if (is.null(dotArgs$libpath) & is.null(dotArgs$shlib)) {
      if(length(dotArgs$useINLAprecomp)>1) {
        warning("length(useINLAprecomp)>1, first taken!")
        dotArgs$useINLAprecomp <- dotArgs$useINLAprecomp[1]
      }
      if(INLAversion>="26.08.22") {
        shlib <-
          cgeneric_shlib_path(
            package = "INLAspacetime",
            useINLAprecomp = FALSE,
            debug = dotArgs$debug
          )
      } else {
        shlib <- cgeneric_shlib_path(
          package = "INLAspacetime",
          useINLAprecomp = dotArgs$useINLAprecomp,
          debug = dotArgs$debug
        )
      }
      if (dotArgs$useINLAprecomp)
        hasverbose <- (INLAversion<="25.02.10") ## to work with old C versions
    } else {
      if(is.null(dotArgs$shlib) & (!is.null(dotArgs$libpath)))
        dotArgs$shlib <- dotArgs$libpath
      shlib <- dotArgs$shlib
      hasverbose <- FALSE
    }

    if(inherits(mesh, "fm_mesh_2d") | inherits(mesh, "inla.mesh")){
      n_mesh <- 1L
      mesh_n <- mesh$n
      bfem <- mesh2fem.barrier(mesh, barrier.triangles)
    } else {
      mesh_n <- sapply(mesh[[1]], function(x) x$n)
      n_mesh <- length(mesh_n)
      bfem <- collect2fem.barrier(mesh, barrier.triangles)
    }
    n <- nrow(bfem$I)
    ndom <- length(bfem$D)

    if(!is.list(barrier.triangles)) {
      barrier.triangles <- list(barrier.triangles)
    }
    if(length(range.fraction) == 1) {
      range.fraction <- rep(range.fraction, ndom-1)
    } else {
      stopifnot(length(range.fraction)==(ndom-1))
    }

    Imat <- bfem$I
    Dmat <- bfem$D[[1]]
    CC <- bfem$C[[1]]
    for(o in 2:ndom) {
      CC <- CC + bfem$C[[o]] * (range.fraction[o-1]^2)
      Dmat <- Dmat +  bfem$D[[o]] * (range.fraction[o-1]^2)
    }
    iC <- Diagonal(n, 1 / CC)
    if(dotArgs$debug) {
      print(utils::str(list(Imat = Imat, Dmat = Dmat, CC = CC, iC = iC)))
    }

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

    if(dotArgs$debug) {
      print(utils::str(lmats))
    }

    args0 <- list(
      model = "inla_cgeneric_barrier",
      shlib = shlib,
      n = as.integer(n),
      debug = as.integer(dotArgs$debug)
    )
    if(hasverbose) { ## to work with old C versions
      args0$verbose <- as.integer(0)
    }
    if(dotArgs$useINLAprecomp && (INLAversion<="25.02.10")) {
      args0$prs <- prior.range
    } else {
      args0$prange <- prior.range
    }
    args0$psigma <- prior.sigma

    the_model <- do.call(
      "cgenericBuilder",
      c(args0,
        list(
        ii = lmats$graph@i,
        jj = lmats$graph@j,
        xx = t(lmats$xx)
      )
      )
    )
    if (constr) {
      At <- matrix(0, n, n_mesh) ## transposed
      idk <- split(1:n, factor(rep(1:n_mesh, mesh_n), 1:n_mesh))
      for(k in 1:n_mesh) {
        if(length(idk[[k]])>0) {
          ## within each domain only
          At[idk[[k]], k] <- bfem$C[[1]][idk[[k]]]
        }
      }
      the_model$f$extraconstr <- list(
        A = t(At), e = rep(0, n_mesh)
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
      if(!is.na(packageCheck(
        name = "inlabru",
        minimum_version = "2.13",
        quietly = TRUE
      ))) {
        the_model$mapper <- inlabru::bm_fmesher(mesh)
      }
    }

    return(the_model)
  }
