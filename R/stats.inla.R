#' To retrieve goodness of fit statistics from an inla output object.
#'
#' Exctracts dic, waic, po, cpo, mse, mae, crps and scrps.
#'
#' @param m an inla output object.
#' @param i an index to subset the estimated values.
#' @param y observed to compare against.
#' @param fsummarize the summary function, default `mean`
#' @section Warning:
#'  The po, crps and scrps are computed using the posterior mean and the
#'  sum of the sd of the posterior fitted value plus the inverse of the
#'  likelihood precision posterior mean.
#' @return a vector with the extracted statistics
#' @export
#'
#' @importFrom stats dnorm pnorm complete.cases
stats.inla <- function(m, i=NULL, y, fsummarize=mean) {
    crps.g <- function(y, m, s) {
        md <- y-m
        s/sqrt(pi) - 2*s*dnorm(md/s) + md*(1-2*pnorm(md/s))
    }
    scrps.g <- function(y, m, s) {
        md <- y-m
        -0.5 * log(2*s/sqrt(pi)) -sqrt(pi) *
            (s*dnorm(md/s)-md/2+md*pnorm(md/s))/s
    }
    if(is.null(i))
      i <- 1:length(m$dic$local.dic)
    r <- c(
      dic=fsummarize(m$dic$local.dic[i]),
      waic=fsummarize(m$waic$local.waic[i]),
      po=-fsummarize(dnorm(
        y[i], m$summary.fitted.value$mean[i],
        sqrt(m$summary.fitted.value$sd[i]^2 +
              1/m$summary.hyperpar$mean[1]), log=TRUE)),
      cpo=-fsummarize(log(m$cpo$cpo[i])),
      mse=fsummarize((m$summary.fitted.value$mean[i]-y[i])^2),
      mae=fsummarize(abs(m$summary.fitted.value$mean[i]-y[i])),
      crps=-fsummarize(crps.g(
                 y[i], m$summary.fitted.value$mean[i],
                 sqrt(m$summary.fitted.value$sd[i]^2 +
                 1/m$summary.hyperpar$mean[1]))),
      scrps=-fsummarize(scrps.g(
        y[i],
        m$summary.fitted.val$mean[i],
        sqrt(m$summary.fitted.val$sd[i]^2 +
             1/m$summary.hyperpar$mean[1]))))
    return(r)
}
