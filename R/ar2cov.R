#' Illustrative code to compute the covariance of
#' the second order autoregression (AR2) model.
#'
#' Computes the auto-covariance for given coefficients.
#'
#' @param a1 the first auto-regression coefficient.
#' @param a2 the second auto-regression coefficient.
#' @param k maximum lag for evaluating the auto-correlation.
#' @param useC just a test (to use C code).
#' @section Details:
#'  Let the second order auto-regression model defined as
#'   `x_t + a_1 x_{t-1} + a_2 x_{t-2} = w_t`
#' where `w_t ~ N(0, 1)`.
#' @return the autocorrelation as a vector or matrix, whenever `a1` or `a2` are
#' scalar or vector.
#' @seealso [ar2precision]
#' @export
#' @examples
#' ar2cov(c(-1.7, -1.8), 0.963, k = 5)
#' plot(ar2cov(-1.7, 0.963), type = "o")
ar2cov <- function(a1, a2, k = 30, useC = FALSE) {
  a <- -cbind(a1, a2)
  n <- nrow(a)
  k <- as.integer(k)
  stopifnot(k>0)
  if(useC) {
    r <- matrix(double(k * n), k, n)
    r[1, ] <- a[, 1] / (1 - a[, 2])
    if (k > 1) {
      r[2, ] <- (a[, 1]^2 + a[, 2] - a[, 2]^2) / (1 - a[, 2])
      if(k>2) {
        r <- .C("ar2cov", n, k, a[,1], a[,2],
                r = r,
                PACKAGE = "INLAspacetime")$r
      }
      r <- t(r)
    }
  } else {
    r <- matrix(double(n * k), n, k)
    r[, 1] <- a[, 1] / (1 - a[, 2])
    if (k > 1) {
      r[, 2] <- (a[, 1]^2 + a[, 2] - a[, 2]^2) / (1 - a[, 2])
    }
    if(k>2) {
      for (j in 3:k) {
        r[, j] <- a[, 1] * r[, j - 1] + a[, 2] * r[, j - 2]
      }
    }
  }
  return(drop(r))
}
