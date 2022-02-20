#' Illustrative code for Finite Element matrices of a mesh in 2d domain.
#'
#' Creates a list of matrices object.
#'
#' @param mesh a 2d mesh object.
#' @param order the desired order.
#' @section Warning:
#'  This is just for illustration purpose and one should consider the
#'  efficient function available a the INLA package.
#' @return a list object containing the FE matrices.
#' @export
fem2d <- function(mesh, order=2) {
  heron <- function(x, y) {
    ### function to compute the area of a triangle
    aa <- sqrt((x[2]-x[1])^2 + (y[2]-y[1])^2)
    bb <- sqrt((x[3]-x[2])^2 + (y[3]-y[2])^2)
    cc <- sqrt((x[1]-x[3])^2 + (y[1]-y[3])^2)
    s <- 0.5*(aa+bb+cc)
    sqrt(s*(s-aa)*(s-bb)*(s-cc))
  }
  areapl <- function(xy) {
    n <- nrow(xy)
    abs(0.5*sum(xy[1:n,1]*xy[c(2:n,1),2]-
                  xy[1:n,2]*xy[c(2:n,1),1]))
  }
  stiffness <- function(x, y) {
    d <- rbind(c(x[3]-x[2], x[1]-x[3], x[2]-x[1]),
               c(y[3]-y[2], y[1]-y[3], y[2]-y[1]))
    crossprod(d)/4
  }
  n <- nrow(mesh$loc)
  c0 <- double(n)
  c1aux <- matrix(1, 3, 3) + diag(3)
  ntv <- nrow(mesh$graph$tv)
  g1 <- c1 <- matrix(0, n, n)
  for (j in 1:ntv) {
    it <- mesh$graph$tv[j,]
    h <- heron(mesh$loc[it,1], mesh$loc[it,2])
    c0[it] <- c0[it] + h/3
    c1[it, it] <- c1[it, it] + h*c1aux/12
    g1[it, it] <- g1[it, it] +
      stiffness(mesh$loc[it, 1],
                mesh$loc[it, 2])/h
  }
  res <- list(c0=as(sparseMatrix(i=1:n, j=1:n, x=c0),
                    'dgTMatrix'),
              c1=as(Matrix(c1, sparse=TRUE),
                    'dgTMatrix'),
              g1=as(Matrix(g1, sparse=TRUE),
                    'dgTMatrix'))
  order <- floor(order)
  if (order>1) {
    for (o in 2:order)
      res[[2+o]] <- as(crossprod(
        res$g1/diag(res$c0),
        res[[o+1]]), 'dgTMatrix')
    names(res)[4:(order+2)] <- paste0('g', 2:order)
  }
  res$va <- matrix(c0, ncol=1)
  if (is.null(mesh$SP)) {
    ##mesh$SP <- sp:::SpatialPolygons(
    ##   lapply(1:nrow(mesh$graph$tv), function(j) {
    ##     p <- sp:::Polygon(mesh$loc[mesh$graph$tv[j, ], 1:2])
    ##    sp:::Polygons(list(p), paste(j))
    ##}))
    ##mesh$centroids <- sp:::coordinates(mesh$SP)
    mesh$ta <-  NULL
  } else {
    res$ta <- matrix(sapply(mesh$SP@polygons, function(p)
      p@area), ncol=1)
  }
  return(res)
}
