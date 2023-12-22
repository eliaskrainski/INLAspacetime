#' Illustrative code for Finite Element matrices of a mesh in 2d domain.
#' @rdname mesh2fem
#' @name mesh2fem
#' @aliases mesh2fem
#' @param mesh a 2d mesh object.
#' @param order the desired order.
#' @param barrier.triangles integer index to specify the
#' triangles in the barrier domain
#' @return a list object containing the FE matrices.
#' @export
mesh2fem <- function(mesh, order = 2, barrier.triangles = NULL) {
  if (!is.null(barrier.triangles)) {
    return(mesh2fem.barrier(
      mesh = mesh,
      barrier.triangles = barrier.triangles
    ))
  }
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
    it <- mesh$graph$tv[j, ]
    h <- Heron(mesh$loc[it, 1], mesh$loc[it, 2])
    ta[j] <- h
    c0[it] <- c0[it] + h / 3
    ii[ncg + 1:9] <- it[ii0]
    jj[ncg + 1:9] <- it[jj0]
    c1x[ncg + 1:9] <- h * c1aux / 12
    g1x[ncg + 1:9] <- Stiffness(
      mesh$loc[it, 1],
      mesh$loc[it, 2]
    ) / h
    ncg <- ncg + 9L
  }
  ijx <- which((ii > 0) | (jj > 0))
  res <- list(
    c0 = Matrix::sparseMatrix(i = 1:n, j = 1:n, x = c0, repr = "T"),
    c1 = INLA::inla.as.dgTMatrix(Matrix::sparseMatrix(
      i = ii[ijx], j = jj[ijx], x = c1x[ijx], dims = c(n, n)
    )),
    g1 = INLA::inla.as.dgTMatrix(Matrix::sparseMatrix(
      i = ii[ijx], j = jj[ijx], x = g1x[ijx], dims = c(n, n)
    ))
  )
  order <- floor(order)
  if (order > 1) {
    for (o in 2:order) {
      g1s <- res$g1 %*% Diagonal(n, 1 / c0)
      res[[2 + o]] <- INLA::inla.as.dgTMatrix(g1s %*% res[[o + 1]])
    }
    names(res)[4:(order + 2)] <- paste0("g", 2:order)
    res$va <- matrix(c0, ncol = 1)
    res$ta <- matrix(ta, ncol = 1)
  }
  return(res)
}
#' Illustrative code for Finite Element matrices when some triangles are
#' in a barrier domain.
#' @rdname mesh2fem
#' @name mesh2fem
#' @aliases mesh2fem.barrier
#' @return a list object containing the FE matrices
#' for the barrier problem.
#' @export
mesh2fem.barrier <- function(mesh, barrier.triangles = NULL) {
  if (is.null(barrier.triangles)) {
    warning("No 'barrier.triangles', using 'mesh2fem(mesh, order = 2)'!")
    return(mesh2fem(mesh = mesh, order = 2L))
  }

  stopifnot(fm_manifold(mesh, c("S", "R")))
  Rmanifold <- fm_manifold(mesh, "R") + 0L

  dimension <- fm_manifold_dim(mesh)
  stopifnot(dimension > 0)

  ### safe to use from triangle area from:
  hh <- INLA::inla.mesh.fem(mesh, order = 1)$ta
  if(Rmanifold==0) {
    stopifnot(all.equal(mesh$crs, fmesher:::fm_crs("sphere"))) ### assume sphere
    ll_crs <- sf::st_crs("+proj=longlat  +datum=WGS84")
    po_crs <- sf::st_crs("+proj=moll +x_0=0 +y_0=0 +lat_0=0 +lon_0=180 +units=km")
    ao_crs <- sf::st_crs("+proj=moll +x_0=0 +y_0=0 +lat_0=0 +lon_0=0 +units=km")
    llloc <- fmesher::fm_transform(mesh, crs = ll_crs)$loc
    ploc <- fmesher::fm_transform(mesh, crs = po_crs)$loc
    mesh <- fmesher::fm_transform(mesh, crs = ao_crs)
    ile <- which((llloc[, 1] < (-175)) | (llloc[, 1]>175))
    Rsq <- 6371^2
    Ea <- 4 * pi * Rsq
    hhS <- hh * Ea
  }

  barrier.triangles <- unique(sort(barrier.triangles))
  itv <- list(
    setdiff(
      1:nrow(mesh$graph$tv),
      barrier.triangles
    ),
    barrier.triangles
  )
  ntv <- sapply(itv, length)
  n <- nrow(mesh$loc)

  c1aux <- c(2, 2, 2, 1, 1, 1, 1, 1, 1)
  ii0c <- c(1, 2, 3, 1, 2, 1, 3, 2, 3)
  jj0c <- c(1, 2, 3, 2, 1, 3, 1, 3, 2)
  ii0g <- rep(1:3, each = 3)
  jj0g <- rep(1:3, 3)

  res <- list(
    I = Diagonal(n, x = rep(0.0, n)),
    D = vector("list", 2L),
    C = vector("list", 2L),
    hdim = 0L
  )

  for (o in 1:2) {
    c0 <- double(n)
    jjc <- iic <- integer(ntv[o] * 9L)
    jjg <- iig <- integer(ntv[o] * 9L)
    g1x <- double(ntv[o] * 9L)
    c1x <- double(ntv[o] * 9L)

    ng <- nc <- 0
    for (j in itv[[o]]) {
      it <- mesh$graph$tv[j, ]
      h <- hh[it]
      hS <- hhS[it]
      if(any(it %in% ile)) {
        ##h <- Heron(ploc[it, 1], ploc[it, 2])
        g1x[ng + 1:9] <- Stiffness(
          ploc[it, 1],
          ploc[it, 2]
        ) / hS
      } else {
        ##h <- Heron(mesh$loc[it, 1], mesh$loc[it, 2])
        g1x[ng + 1:9] <- Stiffness(
          mesh$loc[it, 1],
          mesh$loc[it, 2]
        ) / hS
      }

      c0[it] <- c0[it] + h / 3
      iic[nc + 1:9] <- it[ii0c]
      jjc[nc + 1:9] <- it[jj0c]
      c1x[nc + 1:9] <- h * c1aux / 12

      iig[ng + 1:9] <- it[ii0g]
      jjg[ng + 1:9] <- it[jj0g]

      nc <- nc + 9L
      ng <- ng + 9L
    }

    ijxc <- which((iic > 0) | (jjc > 0))
    ijxg <- which((iig > 0) | (jjg > 0))
    res$I <- res$I +
      INLA::inla.as.dgTMatrix(Matrix::sparseMatrix(
        i = iic[ijxc], j = jjc[ijxc], x = c1x[ijxc], dims = c(n, n)
      ))
    res$D[[o]] <- INLA::inla.as.dgTMatrix(Matrix::sparseMatrix(
      i = iig[ijxg], j = jjg[ijxg], x = g1x[ijxg], dims = c(n, n)
    ))
    res$C[[o]] <- c0 ## * 3
    res$hdim <- res$hdim + 1L
  }
  res$I <- INLA::inla.as.dgTMatrix(res$I)
  res$D[[1]] <- INLA::inla.as.dgTMatrix(res$D[[1]])
  res$D[[2]] <- INLA::inla.as.dgTMatrix(res$D[[2]])
  return(res)
}
