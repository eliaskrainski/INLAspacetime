#' To retrieve the ocean map.
#'
#' Download, if not locally available, the ocean shapefile and
#' load it into the R environment.
#'
#' @param local the directory to save the downloaded file.
#' Default is the current working directory.
#' @section Warning:
#'  it would not work if the file is not available
#'  locally and there is no internet connection.
#' @seealso [worldMap()]
#' @return the oceans map
#' @references
#'  \link{https://maps.princeton.edu/catalog/stanford-jb593mf9286}
#' @export
oceans50m <- function(
  local=".",
  zipfile='oceans50m.zip',
  url="https://stacks.stanford.edu/file/druid:jb593mf9286/data.zip") {
  zipfile <- paste0(local, '/', zipfile)
  if (!file.exists(zipfile)) {
    utils::download.file(url, zipfile)
  }
  utils::unzip(zipfile)
  oc <- rgdal::readOGR('.', 'ne_50m_ocean') ##sf::st_read('.', 'ne_50m_ocean')
  system('rm ne_50m_ocean.dbf ne_50m_ocean-fgdc.xml  ne_50m_ocean-iso19110.xml  ne_50m_ocean-iso19139.xml  ne_50m_ocean.prj  ne_50m_ocean.README.html  ne_50m_ocean.shp  ne_50m_ocean.shp.xml  ne_50m_ocean.shx')
  return(oc)
}
