#' To visualize time series over space.
#'
#' Plot each time series over the map centered at the location.
#'
#' @param stdata matrix with the data, each column is a location.
#' @param spatial an object with one of class defined in the sp package.
#' @param group an integer vector indicating to which spatial unit
#' each time series belongs to. Default is NULL and them it is assumed that
#' each time series belongs o each spatial unit.
#' @param nmax.group an integer indicating the maximum number of time series
#' to be plotted over each spatial unit. Default is NULL, so all will be drawn.
#' @param xscale numeric to define a scaling factor in the horizontal direction.
#' @param yscale numeric to define a scaling factor in the vertical direction.
#' @param colour color (may be a vector, one for each time series).
#' Default is NULL and it will generate colors considering the
#' average of each time series.
#' These automatic colors are defined using the `rgb()` function with `alpha=0.5`.
#' It considers the relative rank of each time series mean, `r`.
#' `r` is then used for red,  `1-r` is used for blue and
#' a triangular function, `1-2*|1-r/2|`, is considered for green.
#' That is, time series with mean among the lowest time series averages
#' are shown in blue and those among the highest temperatures are shown in red.
#' The transition from blue to red goes so that
#' the intermediate ones are shown in light green.
#' @param ... further arguments to be passed for the lines function.
#' @details Scaling the times series is needed before drawing it over the map.
#'  The area of the bounding box for the spatial object
#'  divided by the number of locations is the standard scaling factor.
#'  This is further multiplied by the user given \code{xcale} and \code{yscale}.
#' @section Warning:
#'  if there are too many geographical locations, it will not look good
#' @return add lines to an existing plot
#' @export
stlines <- function(stdata, spatial, group=NULL, nmax.group=NULL,
                    xscale=1, yscale=1, colour=NULL, ...) {
  loc <- sp::coordinates(spatial)
  ns <- nrow(loc)
  nt <- nrow(stdata)
  nd <- ncol(stdata)
  if(is.null(group)) {
    stopifnot(nrow(spatial)>=ncol(stdata))
    gspl <- split(1:nd, 1:nd)
  } else {
    stopifnot(ncol(stdata)==length(group))
    gspl <- split(1:nd, factor(group, 1:ns))
  }
  if(is.null(nmax.group)) nmax.group <- nd
  b <- sp::bbox(spatial)
  s0 <- 0.5*sqrt(diff(b[1,])^2 + diff(b[2,])^2)/(ns^0.8)
  z <- scale(stdata, scale=FALSE)
  z <- z/sqrt(mean(z^2, na.rm=TRUE))
  if(is.null(colour)) {
    r <- rank(attr(z, 'scaled:center'))
    u <- (r-0.5)/nd
    colour <- grDevices::rgb(u, 1-2*abs(u-0.5), 1-u, 0.5)
  }
  for(i in 1:ns) {
    xx <- seq(-s0, s0, length=nt)*xscale + loc[i,1]
    nj <- length(gspl[[i]])
    if(nj>0) {
      for(j in 1:min(nj, nmax.group)) {
        yy <- z[, gspl[[i]][j]]*s0*yscale + loc[i,2]
        graphics::lines(xx, yy, col=colour[gspl[[i]][j]], ...)
        if(any(is.na(yy)))
          graphics::points(xx, yy, col=colour[gspl[[i]][j]], ...)
      }
    }
  }
  return(invisible())
}
