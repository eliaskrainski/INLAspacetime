#' The informations on the GHCN stations.
#'
#' To download the stations information data.
#'
#' @param local the directory to save the downloaded file.
#' Default is the current working directory.
#' @param file the filename to download. Default is 'ghcnd-stations.txt'.
#' @section Warning:
#' it takes time to download in case of bad internet connection
#' @seealso \link{https://www.ncei.noaa.gov/pub/data/ghcn/daily/readme.txt}
#' @return it writes the downloaded file locally and return the information
#' as an SpatialPointsDataFrame object
#' @export
ghcnStations <- function(local=".", file='ghcnd-stations.txt') {
  stfl <- 'ghcnd-stations.txt'
  lfile <- paste0(local, '/', file)
  if(!file.exists(lfile)) {
    utils::download.file(paste0('https://www.ncei.noaa.gov/',
                                'pub/data/ghcn/daily/', stfl),
                         lfile)
  }
  ws <- diff(c(0,11,20,30,37,40,71,75,79,85))
  stations <- utils::read.fwf(lfile, ws, comment.char='')
  colnames(stations) <- c('station', 'latitude', 'longitude', 'elevation',
                          'state', 'name', 'gsn', 'hcn/crn', 'wmo')
  sp::coordinates(stations) <- ~ longitude + latitude
  sp::proj4string(stations) <- sp::CRS("+proj=longlat +datum=WGS84")
  return(stations)
}
