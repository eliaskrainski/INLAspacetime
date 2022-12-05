#' Get the world map
#'
#' Consider the world map, fix the wrap around and convert to the
#' Mollweide projection in km distance (default, can be changed).
#'
#' @param crs a CRS object for the target coordinates.
#' Defalt is the Mollweide projection with units in kilometers.
#' @return the world map
#' @references
#'  the world map data is from the maps R package:
#'  \link{https://cran.r-project.org/web/packages/maps/index.html}
#' @examples
#' wmap <- worldMap()
#' sp::plot(wmap)
#' @export
worldMap <- function(crs=sp::CRS("+proj=moll +units=km")) {

    requireNamespace("maps")
    requireNamespace("maptools")
    
### extract the countries world map
    wrld <- maps::map("world", fill=TRUE, plot=FALSE)
    ID <- sapply(strsplit(wrld$names, ":"), "[", 1L)
    
### put into the SpatialPolygons
    wrld_sp <- maptools::map2SpatialPolygons(wrld, ID=ID)

### identify Antarctica and its main polygon
    ii <- which(sapply(wrld_sp@polygons, function(x) x@ID)=="Antarctica")
    jj <- which.max(sapply(wrld_sp@polygons[[ii]]@Polygons, function(p) p@area))
    
### extract the main polygon
    pl0 <- wrld_sp@polygons[[ii]]@Polygons[[jj]]@coords
    n0 <- nrow(pl0)

### add extreme south
    long0 <- c(180-1e-5, seq(180-1e-5, -180+1e-5, length=2*360), -180+1e-5)
    lat0 <- c(pl0[n0-1,2], rep(-90 +1e-5, length(long0)-1))
    pl1 <- rbind(pl0[1:(n0-1), ], cbind(long0, lat0), pl0[n0,,drop=FALSE])

### substitude polygon 
    wrld_sp@polygons[[ii]] <- SpatialPolygons(list(
        Polygons(list(Polygon(cbind(pl1))), '0')))@polygons[[1]]

### fix the wrap around
    sp::proj4string(wrld_sp) <- "+proj=longlat +datum=WGS84"
    wrld_sp1 <- maptools::nowrapSpatialPolygons(wrld_sp, offset=180)

### projection 
    result <- sp::spTransform(wrld_sp1, crs)
    return(result)

}
