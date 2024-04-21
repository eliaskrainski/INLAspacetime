#' Extracts the dual of a mesh object.
#' @aliases mesh.dual
#' @param mesh a 2d mesh object.
#' @param returnclass if
#' 'list' return a list of polygon coordinates,
#' if "sf" return a 'sf' sfc_multipolygon object,
#' if "sv" return a 'terra', SpatVector object,
#' if "SP" return a 'sp' SpatialPolygons object.
#' @param mc.cores number of threads to be used.
#' @return one of the three in 'returnclass'
#' @export
mesh.dual <- function(mesh,
                      returnclass = c("list", "sf", "sv", "SP"),
                      mc.cores = getOption("mc.cores", 2L)) {
  requireNamespace("parallel")
  if (mesh$manifold != "R2") {
    stop("This only works for R2!")
  }
  crs <- mesh$crs
  returnclass <- match.arg(returnclass)
  ce <- t(sapply(1:nrow(mesh$graph$tv), function(i) {
    colMeans(mesh$loc[mesh$graph$tv[i, ], 1:2])
  }))
  pls <- parallel::mclapply(1:mesh$n, function(i) {
    p <- unique(do.call("rbind", lapply(1:3, function(k) {
      j <- which(mesh$graph$tv[, k] == i)
      if (length(j) > 0) {
        return(rbind(
          ce[j, , drop = FALSE],
          cbind(
            mesh$loc[mesh$graph$tv[j, k], 1] +
              mesh$loc[mesh$graph$tv[j, c(2:3, 1)[k]], 1],
            mesh$loc[mesh$graph$tv[j, k], 2] +
              mesh$loc[mesh$graph$tv[j, c(2:3, 1)[k]], 2]
          ) / 2
        ))
      } else {
        return(ce[j, , drop = FALSE])
      }
    })))
    j1 <- which(mesh$segm$bnd$idx[, 1] == i)
    j2 <- which(mesh$segm$bnd$idx[, 2] == i)
    if ((length(j1) > 0) | (length(j2) > 0)) {
      p <- unique(rbind(
        mesh$loc[i, 1:2], p,
        mesh$loc[mesh$segm$bnd$idx[j1, 1], 1:2] / 2 +
          mesh$loc[mesh$segm$bnd$idx[j1, 2], 1:2] / 2,
        mesh$loc[mesh$segm$bnd$idx[j2, 1], 1:2] / 2 +
          mesh$loc[mesh$segm$bnd$idx[j2, 2], 1:2] / 2
      ))
      yy <- p[, 2] - mean(p[, 2]) / 2 - mesh$loc[i, 2] / 2
      xx <- p[, 1] - mean(p[, 1]) / 2 - mesh$loc[i, 1] / 2
    } else {
      yy <- p[, 2] - mesh$loc[i, 2]
      xx <- p[, 1] - mesh$loc[i, 1]
    }
    p <- p[order(atan2(yy, xx)), ]
    return(switch(
      returnclass,
      list = p,
      sf = sf::st_polygon(list(p[c(1:nrow(p), 1), ])),
      sv = sf::st_polygon(list(p[c(1:nrow(p), 1), ])),
      SP = sp::Polygons(list(sp::Polygon(p)), i)
    ))
  }, mc.cores = mc.cores)
  if(is.na(crs) | is.null(crs))
    crs <- ""
  return(switch(
    returnclass,
    list = pls,
    sf = sf::st_sfc(sf::st_multipolygon(pls),
                    crs = crs),
    sv = terra::vect(sf::st_sfc(sf::st_multipolygon(pls),
                                crs = crs)),
    SP = sp::SpatialPolygons(pls,
                             proj4string = sp::CRS(crs))
  ))
}
