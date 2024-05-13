#' Illustrative code to compute the auto-correlation for an AR2 model.
#'
#' Computes the auto-correlation.
#'
#' @param a1 the first auto-regression coefficient.
#' @param a2 the second auto-regression coefficient.
#' @param k maximum lag for evaluating the auto-correlation.
#' @section Details:
#'  Let the second order auto-regression model defined as
#'   `x_t + a_1 x_{t-1} + a_2 x_{t-2} = w_t`
#' where `w_t ~ N(0, 1)`.
#' @return the autocorrelation as a vector or matrix, whenever `a1` or `a2` are
#' scalar or vector.
#' @seealso [ar2precision]
#' @export
#' @examples
#' plot(ar2cor(-1.7, 0.963), type = "o")
ar2cor <- function(a1, a2, k = 30) {
  a <- -cbind(a1, a2)
  n <- nrow(a)
  k <- as.integer(k)
  stopifnot(k>0)
  r <- matrix(double(n * k), n, k)
  r[, 1] <- a[, 1] / (1 - a[, 2])
  if (k > 1) {
    r[, 2] <- (a[, 1]^2 + a[, 2] - a[, 2]^2) / (1 - a[, 2])
  }
  if (k > 2) {
    if(TRUE) {
      r <- t(.C("ar2cor", n, k, a[,1], a[,2],
                r = t(r),
                PACKAGE = "INLAspacetime")$r)
    } else {
      for (j in 3:k) {
        r[, j] <- a[, 1] * r[, j - 1] + a[, 2] * r[, j - 2]
      }
    }
  }
  return(drop(r))
}
