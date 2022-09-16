#' Define the spacetime model with alpha_t=1, alpha_s=2 and alpha_E=1.
#'
#' This function define the object set for using the 'cgeneric'
#' model tool in INLA.
#'
#' @param tmesh the temporal mesh
#' @param smesh the spatial mesh
#' @param prior.sigma a length two vector for the sigma prior parameters,
#' as U and a so that P(sigma>U)=a. If a=0, then U is taken as, fixed, sigma.
#' @param prior.rs a length two vector for the spatial practical range
#' prior paramters, as U and a so that P(rs<U)=a. If a=0,
#' then U is taken as, fixed, spatial range.
#' @param prior.rt a length two vector for the temporal practical
#' range prior paramters as U and a so that P(rt<U)=a.
#' If a=0 then U is take as, fixed, temporal range. 
#' @details
#' See the paper.
#' @return objects to be used in the cgeneric implementation
#' @export
st121cgeneric <-function(tmesh,
                         smesh,
                         prior.rs,
                         prior.rt,
                         prior.sigma) {

    ### check priors
    stopifnot(length(prior.rs)==2)
    stopifnot(length(prior.rt)==2)
    stopifnot(length(prior.sigma)==2)
    stopifnot(is.numeric(prior.rs))
    stopifnot(is.numeric(prior.rt))
    stopifnot(is.numeric(prior.sigma))

    ### build FE matrices    
    tfe <- inla.mesh.fem(tmesh, 1L)
    sfe <- inla.mesh.fem(smesh, 3L)
    nt <- nrow(tfe$g1)
    stopifnot(nt>1)
    n <- as.integer(nt * nrow(sfe$g1))

### prepare the spacetime matrices
### a) consider 1st order model with 1st basis with mass lumping
    m1 <- Diagonal(nt, c(1.0, rep(0, nt-2), 1.0))
    lmats <- list(  ## g_e^2 (
        dc=kronecker(tfe$c0, sfe$c0),  ## 1       * g_s^6 + g_e^2g_t^2g_s^2
        dg=kronecker(tfe$c0, sfe$g1),  ## 3       * g_s^4 + g_e^2g_t^2
        dg2=kronecker(tfe$c0, sfe$g2), ## 3       * g_s^2
        dg3=kronecker(tfe$c0, sfe$g3), ## 1       * 1
        mc=kronecker(m1, sfe$c0),      ##   g_t   * g_s^4 
        mg1=kronecker(m1, sfe$g1),     ## 2 g_t   * g_s^2 
        mg2=kronecker(m1, sfe$g2),     ##   g_t   * 1 
        hg1=kronecker(tfe$g1, sfe$c0), ##   g_t^2 * g_s^2
        hg2=kronecker(tfe$g1, sfe$g1)) ##   g_t^2 * 1      )

    lmats <- upperPadding(lmats, relative=FALSE)

    llib <- system.file('libs', package='INLAspacetime')

    return(do.call(
        'inla.cgeneric.define',
        list(model="inla_cgeneric_st121_model",
             shlib=paste0(llib, '/cgenericModels.so'),
             n=n, debug=0L,
             prs=prior.rs,
             prt=prior.rt,
             psigma=prior.sigma,
             ii=lmats$graph@i+1L,
             jj=lmats$graph@j+1L,
             xx=t(lmats$xx))))
}

