#' Prepare a matrix or a list of matrices
#' for use in some cgeneric code.
#'
#' Define a graph of the union of the supplied matrices
#' and return the row ordered diagonal plus upper triangle
#' after padding with zeroes each one so that
#' all the returned matrices have the same pattern.
#'
#' @param M a matrix or a list of matrices
#' @param relative logical. If a list of matrices is supplied,
#' it indicates if it is to be returned a relative index
#' and the value for each matrix. See details. 
#' @return If a unique matrix is given, return the
#' upper triangle considering the 'T' representation
#' in the Matrix package. If a list of matrices is given,
#' return a list of two elements: 'graph' and 'xx'.
#' The 'graph' is the union of the graph from each matrix.
#' If relative=FALSE, 'xx' is a matrix with number of column equals
#' the the number of matrices imputed.
#' If relative=TRUE, it is a list of length equal the number
#' of matrices imputed. See details. 
#' @details If relative=FALSE, each colums of 'xx' is the
#' elements of the corresponding matrix after being padded
#' to fill the pattern of the union graph.
#' If relative=TRUE, each element of 'xx' would be a list
#' with a relative index, 'r', for each non-zero elements of
#' each matrix is returned relative to the union graph,
#' and the non-lower elements, 'x', of the corresponding matrix.  
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
#' upperPadding(list(a=A, b=B), relative=TRUE)
upperPadding <- function(M, relative=FALSE) {
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
    graph <- .check(as(Reduce(
      '+',
      lapply(M, function(m) {
        m@x <- m@x*0.0 + 1.0
        return(m)
      })), 'CsparseMatrix'))
    uo <- .uof(graph)
    if(relative) {
        xx <- lapply(M, function(m) {
##            i <- .uof(m)
  ##          r <- pmatch(paste(m@i[i], m@j[i]),
    ##                    paste(graph@i[uo], graph@j[uo]))
            ##      return(list(r=r, x=m@x[i]))
            m <- graph*0 + m
            x <- m@x[uo]
            i <- which(abs(x)>sqrt(.Machine$double.eps))
            return(list(r=i, x=x[i]))
        })
    } else {
        xx <- sapply(M, function(m) {
            m <- graph*0 + m
            return(m@x[uo])
        })
    }
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
