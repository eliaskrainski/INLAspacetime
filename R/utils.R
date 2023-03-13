#' Internal util functions
#' @aliases Heron, Area, Stiffness
#' @param x,y coordinate vectors,
#' @section Warning: Internal functions, not exported.
#' @return the area of a triangle
Heron <- function(x, y) {
  ### function to compute the area of a triangle
  aa <- sqrt((x[2]-x[1])^2 + (y[2]-y[1])^2)
  bb <- sqrt((x[3]-x[2])^2 + (y[3]-y[2])^2)
  cc <- sqrt((x[1]-x[3])^2 + (y[1]-y[3])^2)
  s <- 0.5*(aa+bb+cc)
  sqrt(s*(s-aa)*(s-bb)*(s-cc))
}
#' @return the area of a general polygon
Area <- function(x, y) {
  n <- length(x)
  stopifnot(length(y)==n)
  abs(0.5*sum(x[1:n]*y[c(2:n,1)]-
              y[1:n]*x[c(2:n,1)]))
}
#' @return the stiffness matrix for a triangle
Stiffness <- function(x, y) {
  d <- rbind(c(x[3]-x[2], x[1]-x[3], x[2]-x[1]),
             c(y[3]-y[2], y[1]-y[3], y[2]-y[1]))
  crossprod(d)/4
}
#' @aliases gsConstant
#' @param lgammas the SPDE parameters log(gamma_s, gamma_t, gamma_e)
#' @param smanifold spatial domain.
#' Values could be "S1", "S2", "R1", "R2" and "R3".
#' @return the part of sigma due to spatial constant and gamma_s
gsConstant <- function(lgammas, manifold, alpha) {
  stopifnot(substr(smanifold,1,1) %in% c("R", "S"))
  d <- as.integer(substr(smanifold, 2, 2))
  nu.s <- alpha - d/2
  if(substr(smanifold, 1, 1) == "R") {
    gsCs <- (2*alpha-d)*lgammas[1] +
      lgamma(nu.s) - lgamma(alpha) - (d/2) * log(4*pi)
  }
  if(substr(smanifold, 1, 1) == "S") {
    pi4d <- (4*pi)^(d/2)
    if(lgamma[1]<0.5) {
      gsCs <- (2*alpha-d)*lgammas[1] - 2 * alpha * log(pi4d)
    } else {
      if(d==1) warning("fix this")
      if(lgamma[1]<2) {
        gs2 <- exp(lgammas[1]*2)
        k0 <- 0:60
        gsCs <- log(sum((1 + 2*k0)/(pi4d*(gs2 + k0*(k0+1)))))
      } else {
        gsCs <- 2*(alpha-1)*lgammas[1] - log(alpha-1) - log(pi4d)
      }
    }
  }
  return(gsCs)
}
#' @aliases gammas2params, params2gammas
#' @param alpha.t temporal order of the SPDE
#' @param alpha.s spatial order of the spatial differential operator
#' in the non-separable part.
#' @param alpha.e spatial order of the spatial differential operator
#' in the separable part.
#' @return log(spatial range, temporal range, sigma)
#' @examples
#'  gamma2params(log(c(0, 0, 0)), 1, 2, 1, "R2")
gammas2params <- function(lgammas, alpha.t, alpha.s, alpha.e, smanifold = "R2") {
  alpha <- alpha.e + alpha.s * (alpha.t - 0.5)
  stopifnot(substr(smanifold,1,1) %in% c("R", "S"))
  d <- as.integer(substr(smanifold, 2, 2))
  nu.s <- alpha - d/2
  nu.t <- min(alpha_t - 0.5, nu.s / alpha.s)
  lct <- lgamma(alpha_t-0.5) - lgamma(alpha.t) -0.5 * log(4*pi)
  lgsCs <- gsConstant(lgammas, manifold, alpha)
  return(c(lrs = 0.5 * log(8 * nu.s) - lgammas[1],
           lrt = 0.5 * log(8 * nu.t) - alpha.s * lgammas[1] + lgammas[2],
           lsigma = lct + lgsCs + lgammas[2] - 2*lgammas[3]))
}
#' @param lparams log(spatial range, temporal range, sigma)
#' @return log(gamma_s, gamma_t, gamma_e)
#' @examples
#'  params2gammas(log(c(1, 1, 1)), 1, 2, 1, "R2")
params2gammas <- function(lparams, alpha.t, alpha.s, alpha.e, smanifold = "R2") {
  alpha <- alpha.e + alpha.s * (alpha.t - 0.5)
  stopifnot(substr(smanifold,1,1) %in% c("R", "S"))
  d <- as.integer(substr(smanifold, 2, 2))
  nu.s <- alpha - d/2
  nu.t <- min(alpha_t - 0.5, nu.s / alpha.s)
  lct <- lgamma(alpha_t - 0.5) - lgamma(alpha.t) -0.5 * log(4*pi)
  lcRd <- lgamma(alpha_s - d/2) - lgamma(alpha.s) -(d/2) * log(4*pi)
  lg <- c(rs = 0.5 * log(8 * nu.s),
          rt = -0.5 * log(8 * nu.t),
          sigma = 0.5 * (lct + lcRd))
  lg[1] <- lg[1] - lparams[1]
  lg[2] <- lg[2] + alpha.s * lg[1] + lparams[2]
  if(substr(smanifold,1,1)=="s") {
    if(lg[2]<log(.5)) {
      warning("S manifold done as R, fix later")
    } else {
      if(lg[2]<2) {
        warning("S manifold done as R, fix later")
      } else {
        warning("S manifold done as R, fix later")
      }
    }
  }
  aaux <- 0.5 * d - a
  lg[3] <- lg[3]  + aaux * lg[1] - 0.5 * lg[2] - lparams[3]
  return(lg)
}
