#' Get the world map
#'
#' Consider the world map, fix the wrap around and convert to the
#' Mollweide projection in km distance (default, can be changed).
#'
#' @param crs a CRS object for the target coordinates.
#' Defalt is the Mollweide projection with units in kilometers.
#' @param scale see the help from [rnaturalearth::ne_countries()]
#' @param returnclass defined by the [sp] if returnclass = "sp" or
#' by the [sf] package if returnclass = "sf".
#' @return the world map
#' @references
#'  the world map data is from the maps R package.
#' @export
worldMap <- function(crs = sp::CRS("+proj=moll +units=km"),
                     scale = "medium",
                     returnclass = c("sp", "sf")) {
  requireNamespace("sf")
  requireNamespace("rnaturalearth")
  returnclass <- match.arg(returnclass)
  world.ll <- rnaturalearth::ne_countries(
    scale = scale,
    returnclass = returnclass)
  if(returnclass == "sf") {
    result <- sf::st_transform(
      x = world.ll,
      crs = crs)
  } else {
    if(class(crs) == "character")
      crs <- sp::CRS(crs)
    result <- sp::spTransform(
      x = world.ll,
      CRSobj = crs
    )
  }
  return(result)
}
