#' The informations on the GHCN stations.
#'
#' To download the stations information data.
#'
#' @param year an integer to specify the year
#' @param local the directory to save the downloaded file.
#' Default is the current working directory.
#' @section Warning:
#' it takes time to download in case of bad internet connection
#' @seealso \link{https://www.ncei.noaa.gov/pub/data/ghcn/daily/readme.txt}
#' @return it writes the downloaded file locally and return the information
#' as an SpatialPointsDataFrame object
#' @export
ghcnStations <- function(local=".", file='ghcdn-stations.txt') {
  stfl <- 'ghcnd-stations.txt'
  lfile <- paste0(local, '/', file)
  if(!file.exists(lfile)) {
    download.file(paste0('https://www.ncei.noaa.gov/',
                         'pub/data/ghcn/daily/', stfl),
                  lfile)
  }
  ws <- diff(c(0,11,20,30,37,40,71,75,79,85))
  stations <- read.fwf(lfile, ws, comment.char='')
  colnames(stations) <- c('station', 'latitude', 'longitude', 'elevation',
                          'state', 'name', 'gsn', 'hcn/crn', 'wmo')
  coordinates(stations) <- ~ longitude + latitude
  stations@proj4string <- CRS("+proj=longlat +datum=WGS84")
  return(stations)
}
