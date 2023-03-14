#' Illustrative code for Finite Element matrices of a mesh in 2d domain.
#' @aliases mesh2fem
#' @param mesh a 2d mesh object.
#' @param order the desired order.
#' @param barrierTriangles integer index to specify the
#' triangles in the barrier domain
#' @return a list object containing the FE matrices.
#' @export
mesh2fem <- function(mesh, order=2, barrierTriangles = NULL) {
  if(!is.null(barrierTriangles))
    return(mesh2fem.barrier(mesh, order, barrierTriangles))
  n <- nrow(mesh$loc)
  ntv <- nrow(mesh$graph$tv)
  ta <- rep(0, ntv)
  c1aux <- matrix(1, 3, 3) + diag(3)
  ii0 <- rep(1:3, each = 3)
  jj0 <- rep(1:3, 3)
  c0 <- double(n)
  jj <- ii <- integer(n * 9L)
  g1x <- c1x <- double(n * 9L)
  ncg <- 0
  for (j in 1:ntv) {
    it <- mesh$graph$tv[j,]
    h <- INLAspacetime:::Heron(mesh$loc[it,1], mesh$loc[it,2])
    ta[j] <- h
    c0[it] <- c0[it] + h/3
    ii[ncg + 1:9] <- it[ii0]
    jj[ncg + 1:9] <- it[jj0]
    c1x[ncg + 1:9] <- h * c1aux/12
    g1x[ncg + 1:9] <- INLAspacetime:::Stiffness(
      mesh$loc[it, 1],
      mesh$loc[it, 2])/h
    ncg <- ncg + 9L
  }
  ijx <- which((ii>0) | (jj>0))
  res <- list(
    c0=Matrix::sparseMatrix(i=1:n, j=1:n, x=c0, repr = "T"),
    c1=INLA::inla.as.dgTMatrix(Matrix::sparseMatrix(
      i = ii[ijx], j = jj[ijx], x = c1x[ijx], dims = c(n, n))),
    g1=INLA::inla.as.dgTMatrix(Matrix::sparseMatrix(
      i = ii[ijx], j = jj[ijx], x = g1x[ijx], dims = c(n, n)))
  )
  order <- floor(order)
  if(order>1) {
    for (o in 2:order) {
      g1s <- res$g1 %*% Diagonal(n, 1/c0)
      res[[2+o]] <- INLA::inla.as.dgTMatrix(g1s %*% res[[o+1]])
    }
    names(res)[4:(order+2)] <- paste0('g', 2:order)
    res$va <- matrix(c0, ncol = 1)
    res$ta <- matrix(ta, ncol = 1)
  }
  return(res)
}
#' @aliases mesh2fem.barrier
#' @return a list object containing the FE matrices
#' for the barrier problem.
#' @export
mesh2fem.barrier <- function(mesh, order=2, barrierTriangles = NULL) {
  if(is.null(barrierTriangles))
    return(mesh2fem(mesh, order))
  barrierTriangles <- unique(sort(barrierTriangles))
  itv <- list(setdiff(
    1:nrow(mesh$graph$tv),
    barrierTriangles),
    barrierTriangles)
  ntv <- sapply(itv, length)
  n <- nrow(mesh$loc)

  c1aux <- c(2,2,2, 1,1,1, 1,1,1)
  ii0c <- c(1,2,3, 1,2,1,3, 2,3)
  jj0c <- c(1,2,3, 2,1,3,1, 3,2)
  ii0g <- rep(1:3, each = 3)
  jj0g <- rep(1:3, 3)

  res <- list(
    I = Diagonal(n, x = rep(0.0, n)),
    D = vector("list", 2L),
    C = vector("list", 2L),
    hdim = 0L)

  for(o in 1:2) {
    c0 <- double(n)
    jjc <- iic <- integer(ntv[o] * 9L)
    jjg <- iig <- integer(ntv[o] * 9L)
    g1x <- double(ntv[o] * 9L)
    c1x <- double(ntv[o] * 9L)

    ng <- nc <- 0
    for (j in itv[[o]]) {
      it <- mesh$graph$tv[j,]
      h <- INLAspacetime:::Heron(mesh$loc[it,1], mesh$loc[it,2])

      c0[it] <- c0[it] + h/3
      iic[nc + 1:9] <- it[ii0c]
      jjc[nc + 1:9] <- it[jj0c]
      c1x[nc + 1:9] <- h * c1aux/12

      iig[ng + 1:9] <- it[ii0g]
      jjg[ng + 1:9] <- it[jj0g]
      g1x[ng + 1:9] <- INLAspacetime:::Stiffness(
        mesh$loc[it, 1],
        mesh$loc[it, 2])/h

      nc <- nc + 9L
      ng <- ng + 9L
    }

    ijxc <- which((iic>0) | (jjc>0))
    ijxg <- which((iig>0) | (jjg>0))
    res$I <- res$I +
      INLA::inla.as.dgTMatrix(Matrix::sparseMatrix(
        i = iic[ijxc], j = jjc[ijxc], x = c1x[ijxc], dims = c(n, n)))
    res$D[[o]] <- INLA::inla.as.dgTMatrix(Matrix::sparseMatrix(
      i = iig[ijxg], j = jjg[ijxg], x = g1x[ijxg], dims = c(n, n)))
    res$C[[o]] <- c0 * 3
    res$hdim <- res$hdim + 1L
  }
  res$I <- INLA::inla.as.dgTMatrix(res$I)
  res$D[[1]] <- INLA::inla.as.dgTMatrix(res$D[[1]])
  res$D[[2]] <- INLA::inla.as.dgTMatrix(res$D[[2]])
  return(res)
}
