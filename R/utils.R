#' Internal util functions
#' @aliases Heron, Area, Stiffness
#' @param x,y coordinate vectors,
#' @section Warning: Internal functions, not exported.
#' @return the area of a triangle
Heron <- function(x, y) {
  ### function to compute the area of a triangle
  aa <- sqrt((x[2]-x[1])^2 + (y[2]-y[1])^2)
  bb <- sqrt((x[3]-x[2])^2 + (y[3]-y[2])^2)
  cc <- sqrt((x[1]-x[3])^2 + (y[1]-y[3])^2)
  s <- 0.5*(aa+bb+cc)
  sqrt(s*(s-aa)*(s-bb)*(s-cc))
}
#' @return the area of a general polygon
Area <- function(x, y) {
  n <- length(x)
  stopifnot(length(y)==n)
  abs(0.5*sum(x[1:n]*y[c(2:n,1)]-
              y[1:n]*x[c(2:n,1)]))
}
#' @return the stiffness matrix for a triangle
Stiffness <- function(x, y) {
  d <- rbind(c(x[3]-x[2], x[1]-x[3], x[2]-x[1]),
             c(y[3]-y[2], y[1]-y[3], y[2]-y[1]))
  crossprod(d)/4
}
#' @aliases gsConstant
#' @param lgammas the SPDE parameters log(gamma_s, gamma_t, gamma_e)
#' @param smanifold spatial domain.
#' Values could be "S1", "S2", "R1", "R2" and "R3".
#' @return the part of sigma due to spatial constant and gamma_s
gsConstant <- function(lgammas, manifold, alpha) {
  stopifnot(substr(smanifold,1,1) %in% c("R", "S"))
  d <- as.integer(substr(smanifold, 2, 2))
  nu.s <- alpha - d/2
  if(substr(smanifold, 1, 1) == "R") {
    gsCs <- (2*alpha-d)*lgammas[1] +
      lgamma(nu.s) - lgamma(alpha) - (d/2) * log(4*pi)
  }
  if(substr(smanifold, 1, 1) == "S") {
    pi4d <- (4*pi)^(d/2)
    if(lgamma[1]<0.5) {
      gsCs <- (2*alpha-d)*lgammas[1] - 2 * alpha * log(pi4d)
    } else {
      if(d==1) warning("fix this")
      if(lgamma[1]<2) {
        gs2 <- exp(lgammas[1]*2)
        k0 <- 0:60
        gsCs <- log(sum((1 + 2*k0)/(pi4d*(gs2 + k0*(k0+1)))))
      } else {
        gsCs <- 2*(alpha-1)*lgammas[1] - log(alpha-1) - log(pi4d)
      }
    }
  }
  return(gsCs)
}
#' @aliases gammas2params, params2gammas
#' @param alpha.t temporal order of the SPDE
#' @param alpha.s spatial order of the spatial differential operator
#' in the non-separable part.
#' @param alpha.e spatial order of the spatial differential operator
#' in the separable part.
#' @return log(spatial range, temporal range, sigma)
#' @examples
#'  gamma2params(log(c(0, 0, 0)), 1, 2, 1, "R2")
gammas2params <- function(lgammas, alpha.t, alpha.s, alpha.e, smanifold = "R2") {
  alpha <- alpha.e + alpha.s * (alpha.t - 0.5)
  stopifnot(substr(smanifold,1,1) %in% c("R", "S"))
  d <- as.integer(substr(smanifold, 2, 2))
  nu.s <- alpha - d/2
  nu.t <- min(alpha_t - 0.5, nu.s / alpha.s)
  lct <- lgamma(alpha_t-0.5) - lgamma(alpha.t) -0.5 * log(4*pi)
  lgsCs <- gsConstant(lgammas, manifold, alpha)
  return(c(lrs = 0.5 * log(8 * nu.s) - lgammas[1],
           lrt = 0.5 * log(8 * nu.t) - alpha.s * lgammas[1] + lgammas[2],
           lsigma = lct + lgsCs + lgammas[2] - 2*lgammas[3]))
}
#' @param lparams log(spatial range, temporal range, sigma)
#' @return log(gamma_s, gamma_t, gamma_e)
#' @examples
#'  params2gammas(log(c(1, 1, 1)), 1, 2, 1, "R2")
params2gammas <- function(lparams, alpha.t, alpha.s, alpha.e, smanifold = "R2") {
  alpha <- alpha.e + alpha.s * (alpha.t - 0.5)
  stopifnot(substr(smanifold,1,1) %in% c("R", "S"))
  d <- as.integer(substr(smanifold, 2, 2))
  nu.s <- alpha - d/2
  nu.t <- min(alpha_t - 0.5, nu.s / alpha.s)
  lct <- lgamma(alpha_t - 0.5) - lgamma(alpha.t) -0.5 * log(4*pi)
  lcRd <- lgamma(alpha_s - d/2) - lgamma(alpha.s) -(d/2) * log(4*pi)
  lg <- c(rs = 0.5 * log(8 * nu.s),
          rt = -0.5 * log(8 * nu.t),
          sigma = 0.5 * (lct + lcRd))
  lg[1] <- lg[1] - lparams[1]
  lg[2] <- lg[2] + alpha.s * lg[1] + lparams[2]
  if(substr(smanifold,1,1)=="s") {
    if(lg[2]<log(.5)) {
      warning("S manifold done as R, fix later")
    } else {
      if(lg[2]<2) {
        warning("S manifold done as R, fix later")
      } else {
        warning("S manifold done as R, fix later")
      }
    }
  }
  aaux <- 0.5 * d - a
  lg[3] <- lg[3]  + aaux * lg[1] - 0.5 * lg[2] - lparams[3]
  return(lg)
}
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
#' \eqn{y_t = \sum_{k=j}^h w_j * (x_{t-j}+x_{t+j})/2}
#'
#' Both `s` and `ss` are used to define outliers if
#'  \eqn{|x_t-m|/s>ff_1} or \eqn{|x_t-y_t|/ss>ff_2}
#'
#' - attr(, 'xs') the smoothed time series \eqn{y_t}
outDetect <- function(x, weights=NULL, ff = c(7,7))
{
  ### detect outliers in a time serifes and compute standard deviation for segments
  ## x: n length time series data
  ## weights: half way weights for the smoothing
  ## ff: factors for detecting outliers (#stdev)
  ## jumps: jump size in the moving window to compute local sd
  n <- length(x)
  if(is.null(weights)) {
    h <- 7
    ## define triangular weights for the time smoothing
    weights <- h:1
  } else {
    stopifnot(all(diff(weights)<=0))
    h <- length(weights)
  }
  weights <- c(rev(weights), 0, weights)
  m <- mean(x, na.rm=TRUE)
  s <- sd(x, na.rm=TRUE)
  x <- x - m
  xx <- c(rep(NA, h), x, rep(NA, h))
  xs <- x*0
  ii <- which(complete.cases(xs))
  if(length(ii)>0) {
    for(i in ii) {### time smoothing
      xw <- weights * xx[-h:h + i + h]
      sw <- sum(weights[complete.cases(xw)])
      xs[i] <- sum(xw, na.rm = TRUE)/sw
    }
  }
  ss <- sd(x - xs, na.rm = TRUE)
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
#' of `st`, where is returned as an attribute
#' ad detailed below.
#' @section Attributes:
#' - attr(, 'st') numeric vector with the
#' standard deviation at each segment of the data.
stdSubs <- function(x, nsub=12, fs=15)
{
  n <- length(x)
  ## compute stdev for each segment of the time series
  st <- sapply(split(x, 0:(n-1)%/%nsub), function(xw) {
    if(mean(is.na(xw))>0.5) return(NA)
    return(sd(xw, na.rm=TRUE))
  })
  st.m <- mean(st, na.rm=TRUE)
  r <- any((st/st.m)>fs, na.rm = TRUE) |
    any((st.m/st)>fs, na.rm = TRUE)
  attr(r, 'st') <- st
  return(r)
}
