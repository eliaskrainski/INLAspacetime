# Computes the Whittle-Matern correlation function.

This computes the correlation function as derived in Matern model, see
Matern (1960) eq. (2.4.7). For nu=1, see Whittle (1954) eq. (68). For
the limiting case of nu=0, see Besag (1981) eq. (14-15).

## Usage

``` r
cWhittleMatern(x, range, nu, kappa = sqrt(8 * nu)/range)
```

## Arguments

- x:

  distance.

- range:

  practical range (our prefered parametrization) given as range = sqrt(8
  \* nu) / kappa, where kappa is the scale parameter in the specialized
  references.

- nu:

  process smoothness parameter.

- kappa:

  scale parameter, commonly considered in the specialized literature.

## Value

the correlation.

## Details

Whittle, P. (1954) On Stationary Processes in the Plane. Biometrika,
Vol. 41, No. 3/4, pp. 434-449. http://www.jstor.org/stable/2332724

Matern, B. (1960) Spatial Variation: Stochastic models and their
application to some problems in forest surveys and other sampling
investigations. PhD Thesis.

Besag, J. (1981) On a System of Two-Dimensional Recurrence Equations.
JRSS-B, Vol. 43 No. 3, pp. 302-309. https://www.jstor.org/stable/2984940

## Examples

``` r
plot(function(x) cWhittleMatern(x, 1, 5),
  bty = "n", las = 1,
  xlab = "Distance", ylab = "Correlation"
)
plot(function(x) cWhittleMatern(x, 1, 1), add = TRUE, lty = 2)
plot(function(x) cWhittleMatern(x, 1, 0.5), add = TRUE, lty = 3)
abline(h = 0.139, lty = 3, col = gray(0.5, 0.5))
```
