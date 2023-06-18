#' To retrieve goodness of fit statistics.
#'
#' Extracts dic, waic and log-cpo from an output returned by the inla function
#' from the INLA package or by the bru function from the inlabru package,
#' and computes log-po, mse, mae, crps and scrps for a given input.
#' A summary is applied considering the user imputed function,
#' which by default is the mean.
#'
#' @param m an inla output object.
#' @param i an index to subset the estimated values.
#' @param y observed to compare against.
#' @param fsummarize the summary function,
#' the default is [base::mean()].
#' @section Details:
#'  It assumes Gaussian posterior predictive distributions!
#'  Considering the defaults, for n observations,
#'  \eqn{y_i, i = 1, 2, ..., n}, we have
#'
#'  . dic \deqn{\sum_i d_i/n}{%
#'  \sum_i d_i /n}
#'  where \eqn{d_i} is the dic computed for observation i.
#'
#'  . waic \deqn{\sum_i w_i/n}{%
#'  \sum_i w_i / n}
#'  where \eqn{w_i} is the waic computed for observation i.
#'
#'  . lcpo \deqn{-\sum_i \log(p_i)/n}{%
#'  -\sum_i \log(cpo_i) / n}
#'  where \eqn{p_i} is the cpo computed for observation i.
#'
#'  For the log-po, crps, and scrps scores it assumes a
#'  Gaussian predictive distribution for each observation
#'  \eqn{y_i} which the following definitions:
#'  \eqn{z_i = (y_i-\mu_i)/\sigma_i},
#'  \eqn{\mu_i} is the posterior mean for the linear predictor,
#'  \eqn{\sigma_i = \sqrt{v_i + 1/\tau_y}},
#'  \eqn{\tau_y} is the observation posterior mean,
#'  \eqn{v_i} is the posterior variance of the
#'  linear predictor for \eqn{y_i}.
#'
#'  Then we consider \eqn{\phi()} the density of a standard
#'  Gaussian variable and \eqn{\psi()} the corresponding
#'  Cumulative Probability Distribution.
#'
#'  . lpo \deqn{-\sum_i \log(\phi(z_i))/n}{%
#'  \sum_i \log(\phi(z_i))/n}
#'
#'  . crps \deqn{\sum_i r_i/n}{%
#'  \sum_i r_i/n}
#'  where \deqn{r_i=\sigma_i/\sqrt{\pi} - 2\sigma_i\phi(z_i) + (y_i-\mu_i)(1-2\psi(z_i))}{
#'  r_i=\sigma_i/\sqrt{\pi} - 2\sigma_i\phi(z_i) + (y_i-\mu_i)(1-2\psi(z_i))}
#'
#'  . scrps \deqn{\sum_i s_i/n}{%
#'  \frac{\sum_i s_i}{n}}
#'  where \deqn{s_i=-\log(2\sigma_i/\sqrt{\pi})/2 -\sqrt{\pi}(\phi(z_i)-\sigma_iz_i/2+z_i\psi(z_i))}{
#'  s_i=-\log(2\sigma_i/\sqrt{\pi})/2 -\sqrt{\pi}(\phi(z_i)-\sigma_iz_i/2+z_i\psi(z_i))}
#'
#' @section Warning:
#'  All the scores are negatively oriented which means
#'  that smaller scores are better.
#' @return A named numeric vector with the extracted statistics.
#' @references
#' Held, L. and Schrödle, B. and Rue, H. (2009).
#' Posterior and Cross-validatory Predictive Checks:
#' A Comparison of MCMC and INLA.
#' Statistical Modelling and Regression Structures pp 91–110.
#' \url{https://link.springer.com/chapter/10.1007/978-3-7908-2413-1_6}.
#'
#' Bolin, D. and Wallin, J. (2022) Local scale invariance
#' and robustness of proper scoring rules. Statistical Science.
#' \doi{10.1214/22-STS864}.
#' @export
#'
#' @importFrom stats dnorm pnorm complete.cases
stats.inla <- function(m, i = NULL, y, fsummarize = mean) {
  crps.g <- function(y, m, s) {
    md <- y - m
    s / sqrt(pi) - 2 * s * dnorm(md / s) + md * (1 - 2 * pnorm(md / s))
  }
  scrps.g <- function(y, m, s) {
    md <- y - m
    -0.5 * log(2 * s / sqrt(pi)) - sqrt(pi) *
      (s * dnorm(md / s) - md / 2 + md * pnorm(md / s)) / s
  }
  if (is.null(i)) {
    i <- 1:length(m$dic$local.dic)
  }
  sigma2.mean <- INLA::inla.emarginal(
    exp,
    m$internal.marginals.hyperpar[[1]])
  r <- c(
    dic = fsummarize(m$dic$local.dic[i]),
    waic = fsummarize(m$waic$local.waic[i]),
    lpo = -fsummarize(dnorm(
      y[i], m$summary.fitted.value$mean[i],
      sqrt(m$summary.fitted.value$sd[i]^2 +
        sigma2.mean),
      log = TRUE
    )),
    lcpo = -fsummarize(log(m$cpo$cpo[i])),
    mse = fsummarize((y[i] - m$summary.fitted.value$mean[i])^2),
    mae = fsummarize(abs(y[i] - m$summary.fitted.value$mean[i])),
    crps = -fsummarize(crps.g(
      y[i], m$summary.fitted.value$mean[i],
      sqrt(m$summary.fitted.value$sd[i]^2 +
        sigma2.mean)
    )),
    scrps = -fsummarize(scrps.g(
      y[i],
      m$summary.fitted.val$mean[i],
      sqrt(m$summary.fitted.val$sd[i]^2 +
        sigma2.mean)
    ))
  )
  return(r)
}
