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
#' @details
#' See the paper.
#' @return objects to be used in the f() formula term in INLA.
#' @export
stModel.define <-
    function(smesh, tmesh, model, control.priors)
{
    stopifnot(model %in% c('102','121','202','220'))

    manifold <- (smesh$manifold=='R2') + 0L

    alphas <- as.integer(strsplit(model, '')[[1]])
    nu.t <- alphas[1]-1/2
    alpha <- alphas[3] + alphas[2]*nu.t
    nu.s <- alpha-1
    
    cc <- c(0.5*log(8*nu.s), -0.5*log(8*nu.t),
            0.5*(lgamma(nu.t) - lgamma(alphas[1]) -1.5*log(4*pi)))
    if(manifold) {
        cc[3] <- 0.5*(lgamma(nu.t) + lgamma(nu.s) -
                       lgamma(alphas[1]) -lgamma(alpha) -1.5*log(4*pi))
    }
    
    mm <- stModel.matrices(smesh, tmesh, model)
    n <- smesh$n * tmesh$n
    nm <- ncol(mm$TT)
        stopifnot(nm == length(mm$bb))
        jmm <- pmatch(paste0('M', 1:nm), names(mm))
        stopifnot(length(jmm[complete.cases(jmm)]) == nm)

        llib <- system.file("libs", package = "INLAspacetime")
        lmats <- upperPadding(mm[jmm], relative = FALSE)
        stopifnot(n == nrow(lmats$graph))

        return(do.call(
            "inla.cgeneric.define",
            list(model = "inla_cgeneric_sstspde",
                 shlib = paste0(llib, "/cgenericModels"),
                 n = n, debug = 0L,
                 ii = lmats$graph@i,
                 jj = lmats$graph@j,
                 aaa = alphas,
                 manifold = as.integer(manifold),
                 nm = as.integer(nm),
                 cc = as.double(cc),
                 bb = mm$bb,
                 prs = control.priors$prs,
                 prt = control.priors$prt,
                 psigma = control.priors$psigma,
                 tt = t(mm$TT),
                 xx = t(lmats$xx))))
}

