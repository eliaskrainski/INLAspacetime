#' Helper functions to retrieve the world map,
#' a world polygon, and create grid centers.
#'
#' Retrieve the map of the countries
# ' from the 'rnaturalearth' package.
#'
#' @param crs a string with the projection.
#' Default is the Mollweide projection
#' with units in kilometers.
#' @param scale The scale of map to return.
#' Please see the help of 'ne_countries' function
#' from the 'rnaturalearth' package.
#' @param returnclass A string determining the class
#' of the spatial object to return.
#' Please see the help of 'ne_countries' function
#' from the 'rnaturalearth' package.
#' @references
#'  The land and ocean maps are obtained with
#'  the 'rnaturalearth' package.
#' @export
worldMap <- function(crs = "+proj=moll +units=km",
                     scale = "medium",
                     returnclass = c("sf", "sv")) {
  returnclass <- match.arg(returnclass)
  world.ll <- rnaturalearth::ne_countries(
    scale = scale,
    returnclass = returnclass)
  if(returnclass == "sf") {
    result <- sf::st_transform(
      x = world.ll,
      crs = crs)
  } else {
    result <- terra::project(
      x = world.ll,
      y = crs
    )
  }
  return(result)
}
#' Function to define the boundary Earch polygon
#' in longlat projection for a given resolution.
#' @param resol is the number of subdivisions
#' along the latitude coordinates and half the
#' number of subdivisions along the longitude coordinates.
#' @param crs a string with the projection.
#' Default is the Mollweide projection with units in kilometers.
#' @return a 'st_sfc' object with the Earth polygon.
#' @export
Earth_poly <- function(resol = 300,
                       crs = "+proj=moll +units=km") {
  resol <- ceiling(resol)
  stopifnot(resol>1)
  crs0 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  a <- 180 -1e-7
  b <- 90 -1e-7
  x <- seq(-1, 1, length.out = resol * 2 + 1) * a
  y <- seq(-1, 1, length.out = resol + 1) * b
  ret <- sf::st_sfc(sf::st_multipolygon(
    list(sf::st_polygon(list(cbind(
      long = c(x, rep(a, resol-1),
               rev(x), rep(-a, resol)),
      lat = c(rep(b, 2 * resol + 1), rev(y[2:resol]),
              rep(-b, 2 * resol + 1), y[-1])
#      long = c(rev(x), rep(-a, resol-1),
 #              x, rep(a, resol)),
  #    lat = c(rep(-b, 2 * resol + 1), y[2:resol],
   #           rep(b, 2 * resol + 1), rev(y)[-1])
    ))))),
    crs = crs0)
  if(is.null(crs)) return(ret)
  return(sf::st_transform(
    x = ret,
    crs = crs))
}
#' Define a regular grid in 'Mollweide' projection,
#' with units in kilometers.
#' @param size the (in kilometers) of the grid cells.
#' @param domain if provided it should be an `sf` or `sfc` object.
#' In this case, the grid cells with centers falling inside
#' will be retrieved.
#' @return a 'sf' points object with the centers of a
#' grid set within Earth (and the supplied domain)
#' @export
world_grid <- function(size = 50,
                       domain) {
  requireNamespace("sf")
  crs <- "+proj=moll +units=km"
  if(!missing(domain)) {
    if (!inherits(domain, c("sf", "sfc"))) {
      stop("'domain' must be a 'sf' or a 'sfc' object!")
    }
    dcrs <- sf::st_crs(domain)
  }
  Epoly <- sf::st_transform(
    x = Earth_poly(resol = 1000),
    crs = crs)
  bb <- sf::st_bbox(Epoly)
  grid_0 <- sf::st_as_sf(
    expand.grid(x = seq(from = bb[1],
                        to = bb[3],
                        by = size),
                y = seq(from = bb[2],
                        to = bb[4],
                        by = size)),
    coords = 1:2,
    crs = crs)
  ig_sel <- which(sapply(sf::st_intersects(
    grid_0, Epoly), length) > 0)
  if (!missing(domain)) {
    if(sf::st_crs(Epoly) != sf::st_crs(domain)) {
      domain <- sf::st_transform(
        x = domain,
        crs = crs)
    }
    id_sel <- which(sapply(sf::st_intersects(
      grid_0, domain), length) > 0)
    if(length(id_sel)==0) {
      stop("There is no grid centers inside 'domain'!")
    }
    ig_sel <- intersect(ig_sel, id_sel)
  }
  return(grid_0[ig_sel, ])
}
