#' @title Download daily data
#' @description This function helps to download data from NCDC/NOAA.
#' It first checks if the file exists locally.
#' @param year an integer to specify the year
#' @param local the directory to save the downloaded file.
#' Default is the current working directory.
#' @section Warning:
#' it takes time to download in case of bad internet connection
#' @seealso \code{\link{gzVariableSelect}}
#' @return it writes the downloaded file locally
#' @export
ghcnDownload <- function(year, local=".") {
  if(length(year)>1) {
    year <- paste(year[1])
    warning('Only the first element of "year" is taken.')
  }
  if(nchar(year)!=4) {
    stop('"year" must be given in four character format.')
  }
  ufl <- paste0(year, '.csv.gz')
  lfl <- paste0(local, '/', ufl)
  if (file.exists(lfl)) {
    warning(paste('The file', lfl, 'already exists. Please check!'))
  } else {
      download.file(
        paste0('ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/',
               'daily/by_year/', ufl), lfl)
  }
  return(lfl)
}
