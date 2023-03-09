#' Define a spacetime model object for the `f()` call.
#'
#' @param smesh a spatial mesh
#' @param tmesh a temporal mesh
#' @param model a three characters string to specify the
#' smoothness alpha (each one as integer) parameters.
#' Currently it considers the `102`, `121`, `202` and `220` models.
#' @param control.priors a named list with parameter priors.
#' E.g. prior.rs, prior.rt and prior.sigma
#' as vectors with length two (U, a) to define the
#' corresponding PC-prior such that
#' P(r_s<U)=a, P(r_t<U)=a or P(sigma>U)=a.
#' If a=0 then U is taken to be the fixed value of the parameter.
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
stModel.define <-
    function(smesh, tmesh, model, control.priors,
             constr = FALSE, debug=FALSE, verbose=FALSE,
             useINLAprecomp=TRUE, libpath=NULL)
{
    stopifnot(model %in% c('102','121','202','220'))

    stopifnot(any(substr(smesh$manifold,1,1)%in%c('S', 'R')))
    Rmanifold <- (substr(smesh$manifold, 1, 1)=='R') + 0L

    dimension <- as.integer(substr(smesh$manifold, 2, 2))
    stopifnot(dimension>0)

    alphas <- as.integer(strsplit(model, "")[[1]])
    alpha <- alphas[3] + alphas[2] * (alphas[1] - 0.5)
    nu.s <- alpha-dimension/2
    nu.t <-  min(alphas[1] - 0.5, nu.s / alphas[2])

    if(verbose) {
	    print(c(alphas = alphas))
	    print(c(alpha=alpha, nu.s=nu.s, nu.t=nu.t))
    }

    log.C.t <- lgamma(alphas[1] - 0.5) - lgamma(alphas[1]) - 0.5 * log(4*pi)
    cc <- c(c1 = 0.5*log(8*nu.s),
	    c2 = -0.5*log(8*nu.t),
	    c3 = NA)
    if(Rmanifold) {
	log.C.Rd <- lgamma(alpha - (dimension * 0.5)) - lgamma(alpha) -
		(dimension * 0.5) * log(4*pi)
        cc[3] <- 0.5 * ( log.C.t + log.C.Rd )
	if(verbose)
		cat("R manifold, cc[3] = ", cc[3], "\n")
    } else {
	log.C.S2.part <- -log(4*pi) ## S1???
	cc[3] <- 0.5*(log.C.t + log.C.S2.part
		     ) ## c3 part for S2 (S1???), to be completed in C
	if(verbose)
		cat("S manifold, cc[3] = ", cc[3], "\n")
    }
	if(verbose) {
		print(c(cc=cc))
	}

    mm <- stModel.matrices(smesh, tmesh, model, constr)
    n <- smesh$n * tmesh$n
    nm <- ncol(mm$TT)
        stopifnot(nm == length(mm$bb))
        jmm <- pmatch(paste0('M', 1:nm), names(mm))
        stopifnot(length(jmm[complete.cases(jmm)]) == nm)

        if(is.null(libpath)) {
            if(useINLAprecomp) {
                libpath <- INLA::inla.external.lib("INLAspacetime")
            } else {
                libpath <- system.file("libs", package = "INLAspacetime")
                if(Sys.info()['sysname']=='Windows') {
                    libpath <- file.path(libpath, 'INLAspacetime.dll')
                } else {
                    libpath <- file.path(libpath, 'INLAspacetime.so')
                }
            }
        }

        lmats <- upperPadding(mm[jmm], relative = FALSE)
        stopifnot(n == nrow(lmats$graph))

        the_model <- do.call(
          "inla.cgeneric.define",
          list(
            model = "inla_cgeneric_sstspde",
            shlib = libpath,
            n = n,
            debug = as.integer(debug),
            verbose = as.integer(verbose),
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
          )
        )
        if(constr)
          the_model$f$extraconstr <- mm$extraconstr
        # Prepend specialised model class identifier, for bru_mapper use:
        class(the_model) <- c("stModel_cgeneric", class(the_model))
        # Add objects needed by bru_get_mapper.stModel_cgeneric:
        # (alternatively, construct the mapper already here, but that would
        # require loading inlabru even when it's not going to be used)
        the_model[["smesh"]] <- smesh
        the_model[["tmesh"]] <- tmesh

        the_model
    }

