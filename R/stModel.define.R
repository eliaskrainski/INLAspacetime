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
#' If a=0, then U is taken as fixed.
#' @details
#' See the paper.
#' @return objects to be used in the cgeneric or rgeneric
#' @export
stModel.define <-
    function(smesh, tmesh, model, control.priors)
{
    stopifnot(model %in% c('102','121','202','220'))
    aaa <- as.integer(strsplit(model, '')[[1]])

    mm <- stModel.matrices(smesh, tmesh, model)
    n <- smesh$n * tmesh$n
    nm <- length(mm$bb)
    stopifnot(ncol(mm$T)==nm)
    jmm <- pmatch(paste0('M', 1:nm), names(mm))
    jmm <- jmm[complete.cases(jmm)]
    stopifnot(length(jmm)==nm)
    cc <- c(mm$c1, mm$c2, mm$c3)
    stopifnot(length(cc)==3)

        llib <- system.file("libs", package = "INLAspacetime")
        lmats <- upperPadding(mm[jmm], relative = FALSE)
        stopifnot(n == nrow(lmats$graph))
        aaa <- as.integer(strsplit(model, '')[[1]])
        return(do.call(
            "inla.cgeneric.define",
            list(model = "inla_cgeneric_sstspde",
                 shlib = paste0(llib, "/cgenericModels.so"),
                 n = n, debug = 0L,
                 ii = lmats$graph@i,
                 jj = lmats$graph@j,
                 aaa = aaa,
                 nm = as.integer(nm),
                 cc = as.double(cc),
                 bb = mm$bb,
                 prs = control.priors$prs,
                 prt = control.priors$prt,
                 psigma = control.priors$psigma,
                 tt = t(mm$TT),
                 xx = t(lmats$xx))))
}

