#' Illustrative code to build the precision matrix for SPDE kind models.
#'
#' Creates a precision matrix as a sparse matrix object.
#' For general code look at the functions in the INLA package.
#'
#' @param kappa the scale parameter.
#' @param fem a list containing the Finite Element matrices.
#' @param alpha the smoothness parameter.
#' @section Warning:
#'  This is just for illustration purpose and one should consider the
#'  efficient function available a the INLA package.
#' @return the precision matrix as a sparse matrix object.
#' @export
spde2precision <- function(kappa, fem, alpha) {
  result <- fem$c0 * kappa^2
  if (alpha == 1) {
    result <- fem$c0 * kappa^2 + fem$g1
  }
  if (alpha == 2) {
    result <- fem$c0 * kappa^4 + 2 * kappa^2 * fem$g1 + fem$g2
  }
  if ((alpha > 1) & (alpha < 2)) {
    lambda <- alpha - floor(alpha)
    b <- matrix(c(1, 0, 0, 1, 1, 0, 1, 2, 1), 3, 3) %*%
      solve(matrix(
        1 / (c(4:2, 3:1, 2:0) + lambda), 3,
        3
      ), 1 / (c(4:2) + lambda - alpha))
    result <- kappa^4 * fem$c0 * b[1] + kappa^2 * fem$g1 * b[2] + fem$g2 * b[3]
  }
  result <- Matrix::forceSymmetric(result)
  return(result)
}
