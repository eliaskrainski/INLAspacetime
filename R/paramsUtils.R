#' Funtions to help converting from/to user/internal parametrization.
#' @rdname paramsUtils
#' @name paramsUtils
#' @aliases gammas2params
#' @aliases params2gammas
#' @param lgammas the SPDE parameters log(gamma.s, gamma.t, gamma.e)
#' @param smanifold spatial domain manifold.
#' @param alpha the resulting spatial order.
#' Values could be "S1", "S2", "R1", "R2" and "R3".
#' @return the part of sigma due to spatial constant and gamma.s
gsConstant <- function(lgammas, alpha, smanifold) {
  stopifnot(substr(smanifold,1,1) %in% c("R", "S"))
  d <- as.integer(substr(smanifold, 2, 2))
  nu.s <- alpha - d/2
  if(substr(smanifold, 1, 1) == "R") {
    gsCs <- (2*alpha-d)*lgammas[1] +
      lgamma(nu.s) - lgamma(alpha) - (d/2) * log(4*pi)
  } else {
    pi.d <- (4*pi)##^(d/2)
    if(lgamma[1]<log(0.5)) {
      gsCs <- -2*alpha*lgammas[1] -log(pi.d)
    } else {
      if(d==1) warning("fix this")
      if(lgamma[1]<2) {
        gs2 <- exp(lgammas[1]*2)
        k0 <- 0:50
        gsCs <- log(sum((1 + 2*k0)/(pi4d*(gs2 + k0*(k0+1)))))
      } else {
        gsCs <- 2*(alpha-1)*lgammas[1] - log(alpha-1) - log(pi4d)
      }
    }
  }
  return(gsCs)
}
#' @rdname paramsUtils
#' @param alpha.t temporal order of the SPDE
#' @param alpha.s spatial order of the spatial differential operator
#' in the non-separable part.
#' @param alpha.e spatial order of the spatial differential operator
#' in the separable part.
#' @return log(spatial range, temporal range, sigma)
#' @export
#' @examples
#'  gammas2params(log(c(0, 0, 0)), 1, 2, 1, "R2")
gammas2params <- function(lgammas, alpha.t, alpha.s, alpha.e, smanifold = "R2") {
  alpha <- alpha.e + alpha.s * (alpha.t - 0.5)
  stopifnot(substr(smanifold,1,1) %in% c("R", "S"))
  d <- as.integer(substr(smanifold, 2, 2))
  lct <- lgamma(alpha.t-0.5) - lgamma(alpha.t) -0.5 * log(4*pi)
  lgsCs <- gsConstant(lgammas, alpha, smanifold)
  return(c(lrs = 0.5 * log(8 * (alpha - d/2)) - lgammas[1],
           lrt = 0.5 * log(8 * (alpha.t-0.5)) - alpha.s * lgammas[1] + lgammas[2],
           lsigma = lct + lgsCs - lgammas[2] - 2*lgammas[3]))
}
#' @rdname paramsUtils
#' @param lparams log(spatial range, temporal range, sigma)
#' @return log(gamma.s, gamma.t, gamma.e)
#' @export
#' @examples
#'  params2gammas(log(c(1, 1, 1)), 1, 2, 1, "R2")
params2gammas <- function(lparams, alpha.t, alpha.s, alpha.e, smanifold = "R2") {
  alpha <- alpha.e + alpha.s * (alpha.t - 0.5)
  stopifnot(substr(smanifold,1,1) %in% c("R", "S"))
  d <- as.integer(substr(smanifold, 2, 2))
  nu.s <- alpha - d/2
  nu.t <- min(alpha.t - 0.5, nu.s / alpha.s)
  lct <- lgamma(alpha.t - 0.5) - lgamma(alpha.t) -0.5 * log(4*pi)
  lcRd <- lgamma(alpha.s - d/2) - lgamma(alpha.s) -(d/2) * log(4*pi)
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
  aaux <- 0.5 * d - alpha
  lg[3] <- lg[3]  + aaux * lg[1] - 0.5 * lg[2] - lparams[3]
  return(lg)
}
