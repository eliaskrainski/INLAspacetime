#' Illustrative code to build the projector matrix for SPDE models.
#'
#' Creates a projector matrix object.
#'
#' @param mesh a 2d mesh object.
#' @param loc a two columns matrix with the locations to project for.
#' @param lattice Unused; feature not supported by this illustration.
#' @param xlim,ylim vector with the boundary limits.
#' @param dims the number of subdivisions over each boundary limits.
#' @section Warning:
#'  This is just for illustration purpose and one should consider the
#'  efficient functions available in the INLA and inlabru packages,
#'  e.g. `inlabru::fm_evaluator`.
#' @return the projector matrix as a list with sparse matrix object at `x$proj$A`..
#' @export
mesh2projector <- function(mesh, loc = NULL, lattice = NULL,
                           xlim = NULL, ylim = NULL, dims = c(100, 100)) {
  proj1f <- function(i, xy) {
    a <- numeric(3)
    for (j in 1:3) {
      a[j] <- Heron(
        c(mesh$loc[mesh$graph$tv[
          i, -j
        ], 1], xy[1]),
        c(mesh$loc[mesh$graph$tv[
          i, -j
        ], 2], xy[2])
      )
    }
    return(a)
  }
  if (is.null(mesh$SP)) {
    mesh$SP <- sp::SpatialPolygons(
      lapply(1:nrow(mesh$graph$tv), function(j) {
        jj <- mesh$graph$tv[j, ]
        p <- sp::Polygon(mesh$loc[c(jj, jj[1]), ])
        sp::Polygons(list(p), paste(j))
      })
    )
  }
  res <- list()
  if (is.null(loc)) {
    if (is.null(xlim)) {
      xlim <- range(mesh$loc[, 1])
    }
    if (is.null(ylim)) {
      ylim <- range(mesh$loc[, 2])
    }
    res$x <- seq(xlim[1], xlim[2], length = dims[1])
    res$y <- seq(ylim[1], ylim[2], length = dims[2])
    res$lattice <- list(dims = dims, x = res$x, y = res$y)
    res$lattice$loc <- as.matrix(expand.grid(res$x, res$y))
    i <- which((res$lattice$loc[, 1] %in% xlim) |
      (res$lattice$loc[, 2] %in% ylim))
    ns <- length(i)
    res$lattice$segm <- list(
      loc = res$lattice$loc[i, ],
      idx = cbind(
        as.integer(1:ns),
        as.integer(c(2:ns, 1))
      ),
      grp = matrix(as.integer(1), ns, 1),
      is.bnd = TRUE, crs = NULL
    )
    attr(res$lattice$segm, "class") <- "inla.mesh.segment"
    res$crs <- NULL
    attr(res$lattice, "class") <- "inla.mesh.lattice"
    res$loc <- NULL
    m <- prod(dims)
    res$proj <- list(
      t = matrix(as.integer(NA), m, 1),
      bary = matrix(0, m, 1)
    )
    o <- sp::over(sp::SpatialPoints(res$lattice$loc), mesh$SP)
    res$proj$bary <- t(sapply(1:nrow(
      res$lattice$loc
    ), function(x) {
      proj1f(o[x], res$lattice$loc[x, ])
    }))
    ok <- rowSums(is.na(res$proj$bary)) == 0
    res$proj$A <- Matrix::sparseMatrix(
      i = rep((1:nrow(res$lattice$loc))[ok], 3),
      j = mesh$graph$tv[o[ok], ],
      x = as.vector(res$proj$bary[ok, ] /
        rowSums(res$proj$bary[ok, ])),
      dims = c(nrow(res$lattice$loc), nrow(mesh$loc)),
      repr = "C"
    )
    res$proj$ok <- ok
  } else {
    res <- list(x = NULL, y = NULL, lattice = NULL, loc = loc)
    o <- sp::over(sp::SpatialPoints(loc), mesh$SP)
    res$proj <- list(t = matrix(o, ncol = 1))
    res$proj$bary <- t(sapply(1:nrow(loc), function(x) {
      proj1f(o[x], res$loc[x, 1:2])
    }))
    ok <- rowSums(is.na(res$proj$bary)) == 0
    res$proj$A <- as(
      Matrix::sparseMatrix(
        i = rep((1:nrow(res$loc))[ok], 3),
        j = mesh$graph$tv[o[ok], ],
        x = as.vector(res$proj$bary[ok, ] /
          rowSums(res$proj$bary[ok, ])),
        dims = c(nrow(loc), nrow(mesh$loc))
      ),
      "dgCMatrix"
    )
    res$proj$ok <- ok
  }
  return(res)
}
