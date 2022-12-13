#' @title Select the specified variable
#' @description Select a specivied variable from a .gz file and
#' save it into a .RData.
#' @param gzfile the daily data file file. E.g. 2021.csv.gz from
#' \url{https://www.ncei.noaa.gov/pub/data/ghcn/daily/by_year/}
#' see references bellow.
#' @param variable a character indicating which variable to select.
#' @param astype function to be used to type the data.
#' Default is the current working directory.
#' @param qflag Extract data matching this quality flag. Default ''.
#' @param verbose logical to indicate if verbose mode is active
#' @section Warning:
#'  it can take time to execute if the data.table package is not available.
#' @references
#' Menne, M., Durre, I., Vose, R., Gleason, B. and Houston, T. (2012)
#' An overview of the global historical climatology network-daily database.
#' Journal of Atmospheric and Oceanic Technology, 897–910.
#' @return array [days, stations, variables] if more than one variable
#' or a matrix [days, stations] if one variable.
#' @export
gzVariableSelect <- function(gzfile, variable, astype=as.integer,
                             qflag='', verbose=TRUE) {
  if(verbose) t0 <- Sys.time()
  if(requireNamespace("data.table", quietly = TRUE)) {
    d <- data.table::fread(gzfile)
  } else {
    if(verbose) warning('"data.table" is not available... it may take a while.')
    d <- utils::read.csv(gzfile)
  }
  if(verbose) {
    cat('readed ', nrow(d), '')
    t1 <- Sys.time()
    cat('t1 =', t1-t0, '\n')
  }
  ii <- which(d$V3 %in% variable)
  ii <- ii[which(d$V6[ii] %in% qflag)]
  d <- d[ii, ]
  gc(reset = TRUE)
  if(verbose) {
    cat('selected', nrow(d), '')
    t2 <- Sys.time()
    cat('t2 =', t2-t1, '\n')
  }
  cnames <- c("day", "station")
  names(d)[2:1] <- cnames
  if(length(variable)==1) {
      w <- tapply(d$V4, d[, cnames], astype)
  } else {
      cnames <- c(cnames, 'variable')
      names(d)[3] <- 'variable'
      w <- tapply(d$V4, d[, cnames], astype)
      w <- w[, , pmatch(variable, dimnames(w)[[3]]), drop = FALSE]
  }
  if(verbose) {
    cat('dim =', dim(w), '')
    t3 <- Sys.time()
    cat('t3 =', t3-t2, '\n')
  }
  return(w)
}
