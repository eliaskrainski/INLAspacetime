#' Precision matrix for a stationary AR2 model.
#'
#' Creates a precision matrix as a sparse matrix object
#' considering the specification stated in Details.
#'
#' @param n integer with the size of the precision matrix.
#' @param a numeric vector with length three with the coefficients.
#' @details
#' Let the second order auto-regression model be defined as
#'
#' \deqn{a_0 x_t + a_1 x_{t-1} + a_2 x_{t-2} = w_t, w_t ~ N(0, 1).}
#'
#' The stationary assumption is to consider the variance of
#' \eqn{x_t} to be the same for all \eqn{t=1,\ldots,n}.
#' This assumption gives the \eqn{n \times n} symmetric
#' precision matrix \eqn{Q} as a sparse matrix with the
#' following non-zero elements:
#'
#'  \eqn{Q_{1,1} = Q_{n,n} = a_0^2}
#'
#'  \eqn{Q_{2,2} = Q_{n-1,n-1} = a_0^2 + a_1^2}
#'
#'  \eqn{Q_{1,2} = Q_{2,1} = Q_{n-1,n} = Q_{n,n-1} = a_0 a_1}
#'
#'  \eqn{Q_{t,t} = q_0 = a_0^2 + a_1^2 + a_2^2, t = 3, 4, ..., n-2}
#'
#'  \eqn{Q_{t,t-1} = Q_{t-1,t} = q_1 = a_1(a_0 + a_2), t = 3, 4, ..., n-1}
#'
#'  \eqn{Q_{t,t-2} = Q_{t-2,t} = q_2 = a_2 a_0, t = 3, 4, ..., n}
#'
#' @returns sparse matrix.
#' @seealso [ar2cov]
#' @export
#' @examples
#' ar2precision(7, c(1, -1.5, 0.9))
ar2precision <- function(n, a) {
  if (n < 1) {
    return(NULL)
  }
  if (n == 1) {
    return(a[1]^2)
  }
  if (n == 2) {
    return(matrix(
      c(a[1]^2, sum(a[1:2]))[c(1, 2, 2, 1)], 2
    ))
  }
  if (n == 3) {
    return(matrix(
      c(
        a[1]^2, a[1] * a[2], a[1] * a[3],
        a[1] * a[2], sum(a[1:2]^2), a[1] * a[2],
        a[1] * a[3], a[1] * a[2], a[1]^2
      ), 3
    ))
  }
  if (n > 3) {
    return(
      INLAtools::Sparse(
        Matrix::sparseMatrix(
          i = c(1:n, 1:(n - 1), 1:(n - 2)),
          j = c(1:n, 2:n, 3:n),
          x = c(
            a[1]^2, a[1]^2 + a[2]^2,
            rep(sum(a^2), max(0, n - 4)),
            a[1]^2 + a[2]^2, a[1]^2,
            a[1] * a[2],
            rep(a[2] * (a[1] + a[3]), n - 3),
            a[1] * a[2],
            rep(a[1] * a[3], n - 2)
          ),
          symmetric = TRUE,
          repr = "T"
    )))
  }
}
