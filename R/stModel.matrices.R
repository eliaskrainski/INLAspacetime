#' Define the spacetime model matrices.
#'
#' This function computes all the matrices needed to build
#' the precision matrix for a non-stationary and non-separable
#' spatio-temporal model, as in ...
#'
#' @param tmesh a mesh object over the time domain.
#' @param smesh a mesh object over the spatial domain.
#' @param model a string identifying the model.
#' So far we have the following models:
#' '102', '121', '202' and '220' models.
#' @details
#' See the paper and the vignettes.
#' @return a list containing
#' 1. a string with the model identification
#' 2. a length three vector with the constants `c1`, `c2` and `c3`
#' 3. the vector `d`
#' 4. the matrix `T`
#' 5. the model matrices `M_1`, ..., `M_m`
#' @export
stModel.matrices <-
    function(smesh, tmesh, model)
{

    stopifnot(inherits(smesh, "inla.mesh"))
    stopifnot(inherits(tmesh, "inla.mesh.1d"))
    stopifnot(nchar(model)==3)
    stopifnot(model %in% c('102', '121', '202', '220'))

    alphas <- as.integer(strsplit(model, "")[[1]])
    nu.t <- alphas[1]-1/2
    alpha <- alphas[3] + alphas[2]*nu.t
    nu.s <- alpha-1
    val <- list(model=model)
    val$c1 = 0.5*log(8*nu.s)
    val$c2 = -0.5*log(8*nu.t)
    val$c3 = 0.5*(lgamma(nu.t) + lgamma(nu.s) -
                  lgamma(alphas[1]) -lgamma(alpha) -1.5*log(4*pi))

    uM <- function(m) { ### extract the upper
        m <- inla.as.dgTMatrix(m)
        i.u <- which(m@i<=m@j)
        m@i <- m@i[i.u]
        m@j <- m@j[i.u]
        m@x <- m@x[i.u]
        return(m)
    }

    if(model=='102') {
        val$bb <- rep(c(1,2,1), 3)
        val$TT <- rbind(rep(2:0, 3), rep(0:2/2,each=3))
    }
    if(model=='121') {
        val$bb <- c(1,3,3,1, 1,2,1, 1,1)
        val$TT <- rbind(c(3:0, 2:0, 1:0),
                        rep(0:2/2, c(4,3,2)))
    }
    if(model=='202') {
        val$bb <- rep(c(1,2,1), 5)
        val$TT <- rbind(rep(2:0, 5), rep(0:4/2,each=3))
    }
    if(model=='220') {
        val$bb <- c(1,4,6,4,1, 1,3,3,1, 1,2,1, 1,1, 1)
        val$TT <- rbind(c(4:0, 3:0, 2:0, 1:0, 0),
                        rep(0:4/2, c(5,4,3,2,1)))
    }

    tfe <- inla.mesh.fem(tmesh, order = 2)
    sfe <- inla.mesh.fem(smesh, order = 4)

    ns <- nrow(sfe$c0)
    nt <- nrow(tfe$c0)
    h <- mean(diff(tmesh$loc))

    if(model%in%c('102', '121')) {

        J0 <- tfe$c0
        J1 <- sparseMatrix(
            i=c(1, nt), j=c(1, nt), x=c(0.5,0.5))
        J2 <- tfe$g1

        val$M1 <- uM(kronecker(J0, sfe$c0))
        val$M2 <- uM(kronecker(J0, sfe$g1))
        val$M3 <- uM(kronecker(J0, sfe$g2))

        if(model=='102') {

            val$M4 <- uM(kronecker(J1, sfe$c0))
            val$M5 <- uM(kronecker(J1, sfe$g1))
            val$M6 <- uM(kronecker(J1, sfe$g2))

            val$M7 <- uM(kronecker(J2, sfe$c0))
            val$M8 <- uM(kronecker(J2, sfe$g1))
            val$M9 <- uM(kronecker(J2, sfe$g2))

        } else {

            val$M4 <- uM(kronecker(J0, sfe$g3))

            val$M5 <- uM(kronecker(J1, sfe$c0))
            val$M6 <- uM(kronecker(J1, sfe$g1))
            val$M7 <- uM(kronecker(J1, sfe$g2))

            val$M8 <- uM(kronecker(J2, sfe$c0))
            val$M9 <- uM(kronecker(J2, sfe$g1))

        }

    }

    if(model%in%c('202', '220')) {

        Jm <- Jmatrices(tmesh)

        val$M1 <- uM(kronecker(Jm$J0, sfe$c0))
        val$M2 <- uM(kronecker(Jm$J0, sfe$g1))
        val$M3 <- uM(kronecker(Jm$J0, sfe$g2))

        if(model == '202') {

            val$M4 <- uM(kronecker(Jm$J1, sfe$c0))
            val$M5 <- uM(kronecker(Jm$J1, sfe$g1))
            val$M6 <- uM(kronecker(Jm$J1, sfe$g2))

            val$M7 <- uM(kronecker(Jm$J2, sfe$c0))
            val$M8 <- uM(kronecker(Jm$J2, sfe$g1))
            val$M9 <- uM(kronecker(Jm$J2, sfe$g2))

            val$M10 <- uM(kronecker(Jm$J3, sfe$c0))
            val$M11 <- uM(kronecker(Jm$J3, sfe$g1))
            val$M12 <- uM(kronecker(Jm$J3, sfe$g2))

            val$M13 <- uM(kronecker(Jm$J4, sfe$c0))
            val$M14 <- uM(kronecker(Jm$J4, sfe$g1))
            val$M15 <- uM(kronecker(Jm$J4, sfe$g2))

        } else {

            val$M4 <- uM(kronecker(Jm$J0, sfe$g3))
            val$M5 <- uM(kronecker(Jm$J0, sfe$g4))

            val$M6 <- uM(kronecker(Jm$J1, sfe$c0))
            val$M7 <- uM(kronecker(Jm$J1, sfe$g1))
            val$M8 <- uM(kronecker(Jm$J1, sfe$g2))
            val$M9 <- uM(kronecker(Jm$J1, sfe$g3))

            val$M10 <- uM(kronecker(Jm$J2, sfe$c0))
            val$M11 <- uM(kronecker(Jm$J2, sfe$g1))
            val$M12 <- uM(kronecker(Jm$J2, sfe$g2))

            val$M13 <- uM(kronecker(Jm$J3, sfe$c0))
            val$M14 <- uM(kronecker(Jm$J3, sfe$g1))

            val$M15 <- uM(kronecker(Jm$J4, sfe$c0))

        }

    }

    if(length(unique(diff(tmesh$loc)))>1) {
        warning('Edge correction for irregular are not OK yet!')
    }

    return(val)

}

#' The 2nd order temporal matrices with boundary correction
#' @export
Jmatrices <- function(tmesh) {

    tfe <- inla.mesh.fem(tmesh, 2L)
    nt <- nrow(tfe$g1)
    h <- mean(diff(tmesh$loc))

    J0 <- tfe$c0
    J0[2,2] <- tfe$c0[2,2]/2
    J0[nt-1,nt-1] <- tfe$c0[nt-1,nt-1]/2

    J1 <- sparseMatrix(
        i=c(1,1,2, nt-1,nt-1,nt),
        j=c(1,2,2, nt-1,nt,nt),
        x=c(5,-1,5, 5,-1,5)/4)

    J2 <- tfe$g1 * 2
    J2[2,2] <- tfe$g1[2,2]
    J2[1,2] <- J2[2,1] <- tfe$g1[1,2]
    J2[nt,nt-1] <- J2[nt-1,nt] <- tfe$g1[nt,nt-1]
    J2[nt-1,nt-1] <- tfe$g1[nt-1,nt-1]

    J3 <- sparseMatrix(
        i=c(1,1,2, nt-1,nt-1,nt),
        j=c(1,2,2, nt-1,nt,nt),
        x=c(2,-2,2, 2,-2,2)/(h^2))

    J4 <- tfe$g2
    J4[1,1] <- tfe$g2[1,1]/3
    J4[nt,nt] <- tfe$g2[nt,nt]/3
    J4[1,2] <- J4[2,1] <- tfe$g2[1,2]/2
    J4[nt-1,nt] <- J4[nt,nt-1] <- tfe$g2[nt-1,nt]/2
    J4[2,2] <- tfe$g2[2,2]*5/7
    J4[nt-1,nt-1] <- tfe$g2[nt-1,nt-1]*5/7

    return(list(J0=J0, J1=J1, J2=J2, J3=J3, J4=J4))
}

