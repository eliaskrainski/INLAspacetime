#' Computes the Whittle-Matern correlation function.
#'
#' This computes the correlation function as derived in Matern model,
#' see Matern (1960) eq. (2.4.7).
#' For nu=1, see Whittle (1954) eq. (68).
#' For the limiting case of nu=0, see Besag (1981) eq. (14-15).
#'
#' @param x distance.
#' @param range practical range (our prefered parametrization) given as
#' range = sqrt(8 * nu) / kappa,
#' where kappa is the scale parameter in the specialized references.
#' @param nu process smoothness parameter.
#' @param kappa scale parameter, commonly considered in the specialized literature.
#' @section Details:
#'
#' Whittle, P. (1954) On Stationary Processes in the Plane.
#' Biometrika, Vol. 41, No. 3/4, pp. 434-449.
#' http://www.jstor.org/stable/2332724
#'
#' Matern, B. (1960) Spatial Variation: Stochastic models and their application to
#' some problems in forest surveys and other sampling investigations. PhD Thesis.
#'
#' Besag, J. (1981) On a System of Two-Dimensional Recurrence Equations.
#' JRSS-B, Vol. 43 No. 3, pp. 302-309. https://www.jstor.org/stable/2984940
#'
#' @return the correlation.
#' @export
#' @examples
#' plot(function(x) cWhittleMatern(x, 1, 5),
#'   bty = "n", las = 1,
#'   xlab = "Distance", ylab = "Correlation"
#' )
#' plot(function(x) cWhittleMatern(x, 1, 1), add = TRUE, lty = 2)
#' plot(function(x) cWhittleMatern(x, 1, 0.5), add = TRUE, lty = 3)
#' abline(h = 0.139, lty = 3, col = gray(0.5, 0.5))
cWhittleMatern <- function(x, range, nu, kappa = sqrt(8 * nu) / range) {
  x0 <- .Machine$double.eps
  if (nu < x0) nu <- x0
  r0 <- exp((1 - nu) * log(2) + nu * log(kappa * x0) -
    lgamma(nu)) * besselK(kappa * x0, nu)
  res <- exp((1 - nu) * log(2) + nu * log(kappa * x) -
    lgamma(nu)) * besselK(kappa * x, nu) / r0
  res[x < x0] <- 1.0
  return(res)
}
