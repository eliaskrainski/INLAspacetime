#' Detect outliers in a time series considering the raw data
#' and a smoothed version of it.
#' @aliases outDetect
#' @param x numeric vector
#' @param weights non-increasing numeric vector used as weights for
#' computing a smoothed vector as a rooling window average.
#' Default is null and then \eqn{w_j} is proportional to j
#' in the equation in the Details below.
#' @param ff numeric length two vector with the factors
#' used to consider how many times the standard deviation
#' one data point is out to be considered as an outlier.
#' @return logical vector indicating if the data is an outlier
#' with attributes as detailed bellow.
#'
#' - attr(, 'm') is the mean of x.
#'
#' - attr(, 's') is the standard devation of x.
#'
#' - attr(, 'ss') is the standard deviation for
#' the smoothed data \eqn{y_t} that is defined as
#'
#' \eqn{y_t = \sum_{k=j}^h w_j * (x_{t-j}+x_{t+j})/2}
#'
#' Both `s` and `ss` are used to define outliers if
#'
#'   \eqn{|x_t-m|/s>ff_1} or \eqn{|x_t-y_t|/ss>ff_2}
#'
#' - attr(, 'xs') the smoothed time series \eqn{y_t}
outDetect <- function(x, weights = NULL, ff = c(7, 7)) {
  ### detect outliers in a time serifes and compute standard deviation for segments
  ## x: n length time series data
  ## weights: half way weights for the smoothing
  ## ff: factors for detecting outliers (#stdev)
  ## jumps: jump size in the moving window to compute local sd
  n <- length(x)
  if (is.null(weights)) {
    h <- 7
    ## define triangular weights for the time smoothing
    weights <- h:1
  } else {
    stopifnot(all(diff(weights) <= 0))
    h <- length(weights)
  }
  weights <- c(rev(weights), 0, weights)
  m <- mean(x, na.rm = TRUE)
  s <- stats::sd(x, na.rm = TRUE)
  x <- x - m
  xx <- c(rep(NA, h), x, rep(NA, h))
  xs <- x * 0
  ii <- which(complete.cases(xs))
  if (length(ii) > 0) {
    for (i in ii) { ### time smoothing
      xw <- weights * xx[-h:h + i + h]
      sw <- sum(weights[complete.cases(xw)])
      xs[i] <- sum(xw, na.rm = TRUE) / sw
    }
  }
  ss <- stats::sd(x - xs, na.rm = TRUE)
  ### check for outliers in the centered and smoothed data
  r <- (abs(x / s) > ff[1]) | (abs((x - xs) / ss) > ff[2])
  attr(r, "m") <- m
  attr(r, "s") <- s
  attr(r, "ss") <- ss
  attr(r, "xs") <- xs
  return(r)
}
#' To check unusual low/high variance segments
#' @aliases stdSubs
#' @param x numeric vector
#' @param nsub number for the segments length
#' @param fs numeric to use for detecting too
#' hight or too low local standard deviations.
#' @return logical indicating if any of the
#' `st` are `fs` times lower/higher the average
#' of `st`, where is returned as an attribute:
#'
#' - attr(, 'st') numeric vector with the
#' standard deviation at each segment of the data.
stdSubs <- function(x, nsub = 12, fs = 15) {
  n <- length(x)
  ## compute stdev for each segment of the time series
  st <- sapply(split(x, 0:(n - 1) %/% nsub), function(xw) {
    if (mean(is.na(xw)) > 0.5) {
      return(NA)
    }
    return(stats::sd(xw, na.rm = TRUE))
  })
  st.m <- mean(st, na.rm = TRUE)
  r <- any((st / st.m) > fs, na.rm = TRUE) |
    any((st.m / st) > fs, na.rm = TRUE)
  attr(r, "st") <- st
  return(r)
}
