#' To retrieve goodness of fit statistics from an inla output object.
#'
#' Exctracts DIC, WAIC and CPO from an inla output and
#' computes PO, MSE, MAE, CRPS and SCRPS for a given input.
#' A summary is applied considering the user imputed function,
#' which by default is the mean. See Details and Warning!
#'
#' @param m an inla output object.
#' @param i an index to subset the estimated values.
#' @param y observed to compare against.
#' @param fsummarize the summary function,
#' the default is [base::mean()].
#' @section Details:
#'  When computing the PO, CRPS, and SCRPS scores it assumes a
#'  Gaussian predictive distribution with mean equal to the
#'  posterior mean and variance equal to the sum of the
#'  posterior fitted value variance and the inverse of the
#'  observation precision posterior mean.
#' @section Warning:
#'  For CPO, PO, CRPS and SCRPS what is returned is the
#'  negative of the applied summary.
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
      DIC=fsummarize(m$dic$local.dic[i]),
      WAIC=fsummarize(m$waic$local.waic[i]),
      nPO=-fsummarize(dnorm(
        y[i], m$summary.fitted.value$mean[i],
        sqrt(m$summary.fitted.value$sd[i]^2 +
              1/m$summary.hyperpar$mean[1]), log=TRUE)),
      nCPO=-fsummarize(log(m$cpo$cpo[i])),
      MSE=fsummarize((m$summary.fitted.value$mean[i]-y[i])^2),
      MAE=fsummarize(abs(m$summary.fitted.value$mean[i]-y[i])),
      nCRPS=-fsummarize(crps.g(
                 y[i], m$summary.fitted.value$mean[i],
                 sqrt(m$summary.fitted.value$sd[i]^2 +
                 1/m$summary.hyperpar$mean[1]))),
      nSCRPS=-fsummarize(scrps.g(
        y[i],
        m$summary.fitted.val$mean[i],
        sqrt(m$summary.fitted.val$sd[i]^2 +
             1/m$summary.hyperpar$mean[1]))))
    return(r)
}
