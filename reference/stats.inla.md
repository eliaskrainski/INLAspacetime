# To retrieve goodness of fit statistics for a specific model class.

Extracts dic, waic and log-cpo from an output returned by the inla
function from the INLA package or by the bru function from the inlabru
package, and computes log-po, mse, mae, crps and scrps for a given
input. A summary is applied considering the user imputed function, which
by default is the mean.

## Usage

``` r
stats.inla(m, i = NULL, y, fsummarize = mean)
```

## Arguments

- m:

  an inla output object.

- i:

  an index to subset the estimated values.

- y:

  observed to compare against.

- fsummarize:

  the summary function, the default is
  [`base::mean()`](https://rdrr.io/r/base/mean.html).

## Value

A named numeric vector with the extracted statistics.

## Details

It assumes Gaussian posterior predictive distributions! Considering the
defaults, for n observations, \\y_i, i = 1, 2, ..., n\\, we have

. dic \$\$\sum_i d_i/n\$\$ where \\d_i\\ is the dic computed for
observation i.

. waic \$\$\sum_i w_i/n\$\$ where \\w_i\\ is the waic computed for
observation i.

. lcpo \$\$-\sum_i \log(p_i)/n\$\$ where \\p_i\\ is the cpo computed for
observation i.

For the log-po, crps, and scrps scores it assumes a Gaussian predictive
distribution for each observation \\y_i\\ which the following
definitions: \\z_i = (y_i-\mu_i)/\sigma_i\\, \\\mu_i\\ is the posterior
mean for the linear predictor, \\\sigma_i = \sqrt{v_i + 1/\tau_y}\\,
\\\tau_y\\ is the observation posterior mean, \\v_i\\ is the posterior
variance of the linear predictor for \\y_i\\.

Then we consider \\\phi()\\ the density of a standard Gaussian variable
and \\\psi()\\ the corresponding Cumulative Probability Distribution.

. lpo \$\$-\sum_i \log(\phi(z_i))/n\$\$

. crps \$\$\sum_i r_i/n\$\$ where \$\$r_i=\sigma_i/\sqrt{\pi} -
2\sigma_i\phi(z_i) + (y_i-\mu_i)(1-2\psi(z_i))\$\$

. scrps \$\$\sum_i s_i/n\$\$ where \$\$s_i=-\log(2\sigma_i/\sqrt{\pi})/2
-\sqrt{\pi}(\phi(z_i)-\sigma_iz_i/2+z_i\psi(z_i))\$\$

## Warning

All the scores are negatively oriented which means that smaller scores
are better.

## References

Held, L. and Schrödle, B. and Rue, H. (2009). Posterior and
Cross-validatory Predictive Checks: A Comparison of MCMC and INLA.
Statistical Modelling and Regression Structures pp 91–110.
<https://link.springer.com/chapter/10.1007/978-3-7908-2413-1_6>.

Bolin, D. and Wallin, J. (2022) Local scale invariance and robustness of
proper scoring rules. Statistical Science.
[doi:10.1214/22-STS864](https://doi.org/10.1214/22-STS864) .
