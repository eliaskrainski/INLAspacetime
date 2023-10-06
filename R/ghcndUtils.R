#' Donload the data files used
#' @aliases donwloadUtilFiles
#' @param data.dir the folder to store the files.
#' @param year the year of the daily weather data.
#' @param force logical indicating if it is to force
#' the download. If FALSE each file will be downloaded
#' if it does not exists locally yet.
#' @return a named character vector with the local file names:
#' daily.data, stations.all, elevation.
#' @export
downloadUtilFiles <- function(data.dir, year = 2022, force = FALSE) {
  ### 1. daily weather data for one year
  ### 2. stations information
  ### 3. ETOPO2 elevation data

  ### base URL
  ghcnd <- "https://www.ncei.noaa.gov/pub/data/ghcn/daily/"

  ### daily weather data for a given year
  dfl <- paste0(year, ".csv.gz")
  loc.dfl <- file.path(data.dir, dfl)
  if (force | (!file.exists(loc.dfl))) {
    utils::download.file(
      url = paste0(ghcnd, "by_year/", dfl),
      destfile = loc.dfl
    )
  }

  ### all the available stations information
  sfl <- "ghcnd-stations.txt"
  loc.sfl <- file.path(data.dir, sfl)
  if (force | (!file.exists(loc.sfl))) {
    utils::download.file(
      url = paste0(ghcnd, sfl),
      destfile = loc.sfl
    )
  }

  ### elevation data
  efl <- "ETOPO2.RData"
  loc.efl <- file.path(data.dir, efl)
  if (force | (!file.exists(loc.efl))) {
    utils::download.file(
      url = paste0(
        "http://leesj.sites.oasis.unc.edu/",
        "FETCH/GRAB/RPACKAGES/", efl
      ),
      destfile = loc.efl
    )
  }

  return(c(
    daily.data = loc.dfl,
    stations.all = loc.sfl,
    elevation = loc.efl
  ))
}
#' Select data from the daily dataset
#' @aliases ghcndSelect
#' @param gzfile the local filename for
#' the daily data file file. E.g. 2023.csv.gz from
#' \url{https://www.ncei.noaa.gov/pub/data/ghcn/daily/by_year/}
#' see references bellow.
#' @param variable string with the variable name(s) to be selected
#' @param station string (vector) with the station(s) to be selected
#' @param qflag a string with quality control flag(s)
#' @param verbose logical indicating if progress is to be printed
#' @param astype function to convert data to a class,
#' default is set to convert the data to integer.
#' @section Details:
#' The default selects TMIN, TAVG and TMAX and
#' return it as integer because the original data is also integer
#' with units in 10 Celcius degrees.
#' @references
#' Menne, M., Durre, I., Vose, R., Gleason, B. and Houston, T. (2012)
#' An overview of the global historical climatology network-daily database.
#' Journal of Atmospheric and Oceanic Technology, 897–910.
#' @section Warning:
#' It can take time to execute if, for example,
#' the data.table package is not available.
#' @return if more than one variable, it returns an array
#' whose dimentions are days, stations, variables.
#' If one variable, then it returns a matrix whose dimentions
#' are days, stations.
#' @export
ghcndSelect <- function(gzfile,
                        variable = c("TMIN", "TAVG", "TMAX"),
                        station = NULL,
                        qflag = "",
                        verbose = TRUE,
                        astype = as.integer) {
  ### this function selects `variable` from the daily dataset
  ### it select data with the given quality control `qfrag`
  ### it can return the selected data in long or wide format

  if (verbose) {
    t0 <- Sys.time()
  }

  ### read the full dataset
  if (requireNamespace("data.table", quietly = TRUE)) {
    d <- data.table::fread(gzfile, data.table = FALSE)
  } else {
    if (verbose) {
      warning("\"data.table\" is not available... it may take a while.")
    }
    d <- utils::read.csv(gzfile)
  }

  if (verbose) {
    cat("readed ", nrow(d), "observations.")
    t1 <- Sys.time()
    print(t1 - t0)
  }

  ### select the variables and qflag
  ii <- which(d$V3 %in% variable)
  if (verbose) {
    cat("Found ", length(ii), "observations on", variable, "")
    t2 <- Sys.time()
    print(t2 - t1)
  }


  ii <- ii[which(d$V6[ii] %in% qflag)]
    d <- d[ii, ]

  if (verbose) {
    cat("Selected ", length(ii), "observations.")
    t3 <- Sys.time()
    print(t3 - t2)
  }

    if(is.null(station)) {
        t4 <- t3
    } else {
        ii <- which(d$V1 %in% station)
        d <- d[ii, ]
        if (verbose) {
            cat("Selected ", length(ii),
                "observations from", length(station),
                "stations.")
            t4 <- Sys.time()
            print(t4 - t3)
        }        
    }

    if(length(ii)==0) return(NULL)

  cnames <- c("day", "station")
  names(d)[2:1] <- cnames
  if (length(variable) == 1) {
    d <- tapply(d[, 4], d[, cnames[2:1]], astype)
  } else {
    cnames <- c(cnames, "variable")
    names(d)[3] <- "variable"
    d <- tapply(d[, 4], d[, cnames[c(2, 1, 3)]], astype)
    d <- d[, , pmatch(variable, dimnames(d)[[3]]), drop = FALSE]
  }

    if (verbose) {
        cat("Wide data dim =", dim(d), "")
        t5 <- Sys.time()
        print(t5 - t4)
    }

  return(d)
}
