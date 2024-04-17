#' Functions to help converting
#' from/to user/internal parametrization.
#' The internal parameters are
#' 'gamma_s, 'gamma_t', 'gamma_E'
#' The user parameters are
#' 'r_s', 'r_t', 'sigma'
#' @rdname paramsUtils
#' @name paramsUtils
#' @aliases gammas2params
#' @aliases params2gammas
#' @param lg.s the logarithm of the SPDE parameter `\gamma_s`
#' @param alpha the resulting spatial order.
#' @param smanifold spatial domain manifold, which could be
#'  "S1", "S2", "R1", "R2" and "R3".
#' @details See equation (23) in the paper.
#' @return the part of `sigma` from the spatial constant and `\gamma_s`.
lgsConstant <- function(lg.s, alpha, smanifold) {
  stopifnot(substr(smanifold, 1, 1) %in% c("R", "S"))
  d <- as.integer(substr(smanifold, 2, 2))
  nu.s <- alpha - d / 2
  if (substr(smanifold, 1, 1) == "R") {
    gsCs <- (d - 2 * alpha) * lg.s +
      lgamma(nu.s) - lgamma(alpha) - (d / 2) * log(4 * pi)
  } else {
    pi.d <- (4 * pi)^(d/2)
    if (lgamma[1] < log(0.5)) {
      gsCs <- -2 * alpha * lg.s - log(pi.d)
    } else {
      if (d == 1) warning("fix this")
      if (lgamma[1] < 2) {
        gs2 <- exp(lg.s * 2)
        k0 <- 0:50
        gsCs <- log(sum((1 + 2 * k0) / (pi.d * (gs2 + k0 * (k0 + 1))^alpha)))
      } else {
        gsCs <- 2 * (alpha - 1) * lg.s - log(alpha - 1) - log(pi.d)
      }
    }
  }
  return(gsCs)
}
#' Convert from user parameters to SPDE parameters
#' @rdname paramsUtils
#' @param lparams log(spatial range, temporal range, sigma)
#' @return log(gamma.s, gamma.t, gamma.e)
#' @details
#' See equations (19), (20) and (21) in the paper.
#' @export
#' @examples
#' params2gammas(log(c(1, 1, 1)), 1, 2, 1, "R2")
params2gammas <- function(lparams, alpha.t, alpha.s, alpha.e, smanifold = "R2") {
  alpha <- alpha.e + alpha.s * (alpha.t - 0.5)
  stopifnot(substr(smanifold, 1, 1) %in% c("R", "S"))
  d <- as.integer(substr(smanifold, 2, 2))
  nu.s <- alpha - d / 2
  nu.t <- min(alpha.t - 0.5, nu.s / alpha.s)
  lg <- c(
    gamma.s = 0.5 * log(8 * nu.s) - lparams[1],
    gamma.t = -0.5 * log(8 * (alpha.t - 0.5)),
    gamma.e = NA
  )
  lg[2] <- lg[2] + alpha.s * lg[1] + lparams[2]
  lct <- lgamma(alpha.t - 0.5) - lgamma(alpha.t) - 0.5 * log(4 * pi)
  lcs <- lgsConstant(lg[1], alpha, smanifold)
  lg[3] <- 0.5 * (lct + lcs - lg[2]) - lparams[3]
  return(lg)
}
#' Convert from SPDE parameters to user parameters
#' @rdname paramsUtils
#' @param lgammas numeric of length 3 with
#' \eqn{log(\gamma_k)}{log(gamma[k])}
#' model parameters. The parameter order is
#'  log(gamma.s, gamma.t, gamma.e)
#' @param alpha.t temporal order of the SPDE
#' @param alpha.s spatial order of the spatial differential operator
#' in the non-separable part.
#' @param alpha.e spatial order of the spatial differential operator
#' in the separable part.
#' @details
#' See equations (19), (20) and (21) in the paper.
#' @return log(spatial range, temporal range, sigma)
#' @export
#' @examples
#' gammas2params(log(c(0, 0, 0)), 1, 2, 1, "R2")
gammas2params <- function(lgammas, alpha.t, alpha.s, alpha.e, smanifold = "R2") {
  alpha <- alpha.e + alpha.s * (alpha.t - 0.5)
  stopifnot(substr(smanifold, 1, 1) %in% c("R", "S"))
  d <- as.integer(substr(smanifold, 2, 2))
  nu.s <- alpha - d / 2
  nu.t <- min(alpha.t - 0.5, nu.s / alpha.s)
  lct <- lgamma(alpha.t - 0.5) - lgamma(alpha.t) - 0.5 * log(4 * pi)
  lgsCs <- lgsConstant(lgammas[1], alpha, smanifold)
  return(c(
    lrs = 0.5 * log(8 * nu.s) - lgammas[1],
    lrt = 0.5 * log(8 * (alpha.t - 0.5)) - alpha.s * lgammas[1] + lgammas[2],
    lsigma = 0.5 * (lct + lgsCs - lgammas[2]) - lgammas[3]
  ))
}
