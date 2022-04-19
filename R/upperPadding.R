#' Prepare a matrix (or a list of matrices)
#' for use in a cgeneric code.
#'
#' The supplied matrix (or matrices) are prepared to
#' follow the required graph and precision pattern.
#'
#' @param matrices a matrix a list of matrices
#' @return If a unique matrix is given, return the
#' upper triangle considering the 'T' representation
#' in the Matrix package.
#'  If a list of matrices is given, return the graph of the
#' corresponding sum of all the matrices and the corresponding
#' values of the matrices.
#' @export
#' @examples
#' A <- sparseMatrix(
#'    i=c(1,1,2,2,3,3,3,4,5,5),
#'    j=c(2,5,1,3,2,4,5,3,1,3), x=1)
#' A
#' upperPadding(A)
#' B <- Diagonal(5, colSums(A))
#' upperPadding(list(a=A, b=B))
upperPadding <- function(M) {
  if(is(M, 'matrix')) {
    return(upperPadding(Matrix(M)))
  }
  if(is(M, 'Matrix')) {
    if((!is(M, 'dgTMatrix')) &
       (!is(M, 'dtTMatrix')))
      M <- as(as(M, 'dgCMatrix'),
              'dgTMatrix')
    upp <- which(M@j>=M@i)
    ord <- order(M@i[upp])
    return(sparseMatrix(
      i=M@i[upp][ord]+1L,
      j=M@j[upp][ord]+1L,
      x=M@x[upp][ord],
      dims=dim(M),
      repr='T'))
  }
  if(is(M, 'list')) {
    graphMx <- list(
        graph=upperPadding(
          Reduce(
            '+',
            lapply(M, function(m) {
                m@x <- m@x*0.0 + 1.0
                return(m)
            }))))
    graphMx$Mx <- mapply(
        function(x) (x + graphMx$graph*0)@x,
        lapply(M, upperPadding))
    return(graphMx)
  }
  warning('nothing done!')
  return(M)
}
