#' Illustrative code to build the precision matrix for an AR2 model.
#'
#' Creates a precision matrix as a sparse matrix object.
#'
#' @param n the size of the model.
#' @param a the length three vector with the coefficients. See details.
#' @section Details:
#'  Le the second order auto-regression model defined as
#'   a_0 x_t + a_1 x_{t-1} + a_2 x_{t-2} = w_t
#' where w_t ~ N(0, 1).
#' @return the precision matrix as a sparse matrix object with edge correction.
#' @seealso \link{ar2cor}
#' @export
ar2q <- function(n, a) {
  if (n<1) return(NULL)
  if (n==1) return(a[1]^2)
  if (n==2) return(matrix(
    c(a[1]^2, sum(a[1:2]))[c(1,2,2,1)], 2))
  if (n==3) return(matrix(
    c(a[1]^2, a[1]*a[2], a[1]*a[3],
      a[1]*a[2], sum(a[1:2]^2), a[1]*a[2],
      a[1]*a[3], a[1]*a[2], a[1]^2), 3))
  if (n>3)
    return(sparseMatrix(
      i=c(1:n, 1:(n-1), 1:(n-2)),
      j=c(1:n, 2:n, 3:n),
      x=c(a[1]^2, a[1]^2+a[2]^2,
          rep(sum(a^2), max(0,n-4)),
          a[1]^2+a[2]^2, a[1]^2,
          a[1]*a[2],
          rep(a[2]*(a[1]+a[3]), n-3),
          a[1]*a[2],
          rep(a[1]*a[3], n-2)),
      symmetric = TRUE))
}
