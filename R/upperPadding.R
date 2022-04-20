#' Prepare a matrix (or a list of matrices)
#' for use in some cgeneric code.
#'
#' Define a graph of the union of the supplied matrices
#' and return the row ordered diagonal plus upper triangle
#' after padding with zeroes each one so that
#' all the returned matrices have the same pattern.
#'
#' @param M a matrix or a list of matrices
#' @return If a unique matrix is given, return the
#' upper triangle considering the 'T' representation
#' in the Matrix package. If a list of matrices is given,
#' return the union of the graph from each matrix and
#' a matrix with the values each one of the matrices.
#' This matrix has number of columns equal the number
#' of matrices given and number of rows equal the number
#' of elements in the upper triangle and diagonal
#' of the union graph. Return M otherwise.
#' @export
#' @examples
#' A <- sparseMatrix(
#'    i=c(1,1,2,3,3,5),
#'    j=c(2,5,3,4,5,5),
#'    x=-(0:5), symmetric=TRUE)
#' A
#' upperPadding(A)
#' B <- Diagonal(nrow(A), -colSums(A))
#' list(a=A, a=B)
#' upperPadding(list(a=A, b=B))
upperPadding <- function(M) {
## TO DO: relative indexing
  .check <- function(m) {
    if(is(m, 'matrix'))
      m <- Matrix(m)
    if(is(m, 'Matrix')) {
      if((!is(m, 'dgTMatrix')) &
         (!is(m, 'dtTMatrix')))
        m <- as(as(as(m, 'CsparseMatrix'),
                   'dgCMatrix'),
                'dgTMatrix')
      return(m)
    }
    stop('The argument is not valid!')
  }
  .uof <- function(m)
    return(intersect(
      order(m@i),
      which(m@j>=m@i)))
  if(is(M, 'list')) {
    M <- lapply(M, .check)
    graph <- Reduce(
      '+',
      lapply(M, function(m) {
        m@x <- m@x*0.0 + 1.0
        return(m)
      }))
    uo <- .uof(.check(as(graph, 'CsparseMatrix')))
    xx <- sapply(M, function(m) {
      m <- graph*0 + m
      return(m@x[uo])
    })
    return(list(
      graph=sparseMatrix(
        i=graph@i[uo]+1L,
        j=graph@j[uo]+1L,
        x=graph@x[uo],
        dims=dim(graph),
        repr='T'),
      xx=xx))
  } else {
    M <- upperPadding(list(M))
    M$graph@x <- drop(M$xx)
    return(M$graph)
  }
}
