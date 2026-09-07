#' Illustrative code for Finite Element matrices of a mesh in 2d domain.
#' @rdname mesh2fem
#' @name mesh2fem
#' @aliases mesh2fem
#' @param mesh a spatial mesh from [fmesher::fm_mesh_2d()],
#' or collection from [fmesher::fm_collect()]
#' @param order the desired order.
#' @param barrier.triangles integer (or list of, or list of list of) indexes
#' to specify the triangles in the barrier domain.
#' For the only one domain (mesh) and only one barrier it may be vector
#' or a length one list.
#' For only one domain an k barriers, it should be a length k list of indexes
#' for the triangles inside each barrier.
#' If `mesh` is a list of domains, a [fmesher::fm_collect()] output,
#' then it should be a list where each element `i` contains a length k
#' list with a index vector for the triangles in the i-th domain.
#' @return a list object containing the FE matrices.
#' @export
mesh2fem <- function(mesh, order = 2, barrier.triangles = NULL) {
  if (!is.null(barrier.triangles)) {
    return(mesh2fem.barrier(
      mesh = mesh,
      barrier.triangles = barrier.triangles
    ))
  }

  if(inherits(mesh, "fm_collect"))
    return(fmesher::fm_fem(mesh))

  stopifnot(fm_manifold(mesh, c("S", "R")))
  Rmanifold <- fm_manifold(mesh, "R") + 0L

  dimension <- fm_manifold_dim(mesh)
  stopifnot(dimension > 1)

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
    if(Rmanifold==0) {
      h <- s2trArea(mesh$loc[it, ])
    } else {
      h <- Heron(mesh$loc[it, 1], mesh$loc[it, 2])
    }
    ta[j] <- h
    c0[it] <- c0[it] + h / 3
    ii[ncg + 1:9] <- it[ii0]
    jj[ncg + 1:9] <- it[jj0]
    c1x[ncg + 1:9] <- h * c1aux / 12
    fa <- flatArea(mesh$loc[it, ])
    g1x[ncg + 1:9] <- Stiffness(mesh$loc[it, ]) / fa
    ncg <- ncg + 9L
  }
  ijx <- which((ii > 0) | (jj > 0))
  res <- list(
    c0 = Matrix::sparseMatrix(
      i = 1:n, j = 1:n, x = c0),
    c1 = Sparse(Matrix::sparseMatrix(
      i = ii[ijx], j = jj[ijx], x = c1x[ijx],
      dims = c(n, n))),
    g1 = Sparse(Matrix::sparseMatrix(
      i = ii[ijx], j = jj[ijx], x = g1x[ijx],
      dims = c(n, n)))
  )
  order <- floor(order)
  if (order > 1) {
    for (o in 2:order) {
      g1s <- res$g1 %*% Diagonal(n, 1 / c0)
      res[[2 + o]] <- Sparse(g1s %*% res[[o + 1]])
    }
    names(res)[4:(order + 2)] <- paste0("g", 2:order)
    res$va <- matrix(c0, ncol = 1)
    res$ta <- matrix(ta, ncol = 1)
  }
  return(res)
}
#' Finite Element matrices with barrier domain.
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
  } else {
    if(is.list(barrier.triangles)) {
      ntv1 <- nrow(mesh$graph$tv)
      for(i in 1:length(barrier.triangles)) {
        if(length(barrier.triangles[[i]])>0) {
          barrier.triangles[[i]] <- unique(sort(barrier.triangles[[i]]))
          stopifnot(all(barrier.triangles[[i]] %in% (1:ntv1)))
        }
      }
      itv <- c(list(
        setdiff(1:ntv1, unlist(barrier.triangles))),
        barrier.triangles)
      ntv <- sapply(itv, length)
    } else {
      barrier.triangles <- unique(sort(barrier.triangles))
      itv <- list(
        setdiff(
          1:nrow(mesh$graph$tv),
          barrier.triangles
        ),
        barrier.triangles
      )
      ntv <- sapply(itv, length)
    }
  }

  stopifnot(fm_manifold(mesh, c("S", "R")))
  Rmanifold <- fm_manifold(mesh, "R") + 0L

  dimension <- fm_manifold_dim(mesh)
  stopifnot(dimension > 1)

  if(TRUE) {
    ### safe to use from triangle area from:
    hh <- fm_fem(mesh, order = 1)$ta
  } else {
    if(Rmanifold==0) {
      R.i <- sqrt(rowSums(mesh$loc^2))
      hh <- sapply(1:nrow(mesh$graph$tv), function(i) {
        it <- mesh$graph$tv[i, ]
        s2trArea(mesh$loc[it, ], R.i[1])
      })
    } else {
      hh <- sapply(1:nrow(mesh$graph$tv), function(i) {
        it <- mesh$graph$tv[i, ]
        Heron(mesh$loc[it, 1], mesh$loc[it, 2])
      })
    }
  }

  n <- nrow(mesh$loc)

  c1aux <- matrix(1, 3, 3) + diag(3)
  ii0 <- rep(1:3, each = 3)
  jj0 <- rep(1:3, 3)

  res <- list(
    I = Diagonal(n, x = rep(0.0, n)),
    D = vector("list", 2L),
    C = vector("list", 2L),
    hdim = 0L
  )

  for (o in 1:length(itv)) {
    c0 <- double(n)
    ii <- integer(ntv[o] * 9L)
    jj <- integer(ntv[o] * 9L)
    g1x <- double(ntv[o] * 9L)
    c1x <- double(ntv[o] * 9L)

    ng <- nc <- 0
    for (j in itv[[o]]) {
      it <- mesh$graph$tv[j, ]
      h <- hh[j]
      fa <- flatArea(mesh$loc[it, ])
      g1x[ng + 1:9] <- Stiffness(mesh$loc[it, ]) / fa
      c0[it] <- c0[it] + h / 3
      c1x[nc + 1:9] <- h * c1aux / 12
      ii[ng + 1:9] <- it[ii0]
      jj[ng + 1:9] <- it[jj0]
      nc <- nc + 9L
      ng <- ng + 9L
    }

    ijx <- which((ii > 0) | (jj > 0))
    res$I <- res$I + Matrix::sparseMatrix(
      i = ii[ijx], j = jj[ijx], x = c1x[ijx],  dims = c(n, n))
    res$D[[o]] <- Sparse(Matrix::sparseMatrix(
      i = ii[ijx], j = jj[ijx], x = g1x[ijx], dims = c(n, n)))
    res$C[[o]] <- c0
    res$hdim <- res$hdim + 1L
  }
  res$I <- Sparse(res$I)
  return(res)
}
#' Finite Element matrices with barrier domain.
#' @rdname mesh2fem
#' @name mesh2fem
#' @aliases mesh2fem.barrier
#' @param collect an output of [fmesher::fm_collect]
#' @return a list object containing the FE matrices
#' for the barrier problem.
#' @export
collect2fem.barrier <- function(collect, barrier.triangles = NULL) {
  if (is.null(barrier.triangles)) {
    warning("No 'barrier.triangles', using fmesher::fm_fem(, order = 2)'!")
    out <- fmesher::fm_fem(collect, order = 2L)
  } else {
    nc <- length(collect[[1]])
    stopifnot(nc == length(barrier.triangles))
    lfe <- vector("list", nc)
    for(k in 1:nc) {
      lfe[[k]] <- mesh2fem.barrier(collect[[1]][[k]], barrier.triangles[[k]])
    }
    out <- list(
      I = Matrix::.bdiag(lapply(lfe, function(x) x$I))
    )
    ndom <- length(lfe[[1]]$D)
    out$D <- lapply(1:ndom, function(d)
        Matrix::.bdiag(lapply(lfe, function(fe) fe$D[[d]])))
    out$C <- lapply(1:ndom, function(d)
      unlist(lapply(lfe, function(fe) fe$C[[d]])))
  }
  return(out)
}
#' Define the index set of mesh centers within barrier domain.
#'
#' @rdname mesh2fem
#' @name mesh2fem
#' @aliases mesh2fem.barrier
#' @param barriers a list where each element is a
#' polygon defining a barrier domain.
#' @return a list with length of one if the mesh is a `fm_mesh_2d` or
#' length equal the length of the `fm_collect` function spaces.
#' Each of it elements is a list with legth equal the number of barrier domains,
#' containing a index vector with the mesh triangle centroids withing the
#' corresponding barrier.
#' @export
barrier_mesh_centroids <- function(mesh, barriers) {
  if(inherits(mesh, "fm_mesh_2d")) {
    nbr <- length(barriers)
    out <- vector("list", nbr)
    for(k in 1:nbr) {
      out[[k]] <- unlist(
        fmesher::fm_contains(
          x = barriers[[k]],
          y = mesh,
          type = "centroid"
        )
      )
    }
  }
  if(inherits(mesh, "fm_collect")) {
    nfs <- length(mesh$fun_spaces)
    out <- vector("list", nfs)
    for(k in 1:nfs) {
      out[[k]] <- barrier_mesh_centroids(
        mesh$fun_spaces[[k]], barriers
      )
    }
  }
  return(out)
}
