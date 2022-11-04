#' Get the world map 
#'
#' Consider the world map, fix the wrap around and convert to the
#' Mollweide projection in km distance (default, can be changed).
#'
#' @param crs the CRS string to transform the coordinates.
#' Defalt is the Mollweide projection with units in kilometers.
#' @seealso \code{\link{oceans50m}}
#' @return the world map
#' @references
#'  the world map data is from the maps R package:
#'  \link{https://cran.r-project.org/web/packages/maps/index.html}
#' @examples
#' wmap <- worldMap()
#' sp::plot(wmap)
#' @export
worldMap <- function(crs=CRS("+proj=moll +units=km")) {
  ### extract the countries world map
  wrld <- maps::map("world", fill=TRUE, plot=FALSE)
  ID <- sapply(strsplit(wrld$names, ":"), "[", 1L)
  ### put into the SpatialPolygons
  wrld_sp <- maptools::map2SpatialPolygons(wrld, ID=ID)
  proj4string(wrld_sp) <- "+proj=longlat +datum=WGS84"
  ### fix the wrap around
  wrld_sp1 <- maptools::nowrapSpatialPolygons(wrld_sp, offset=180)
  ### convert map from longlat to mollweide
  result <- spTransform(wrld_sp1, crs)
  return(result)
}
