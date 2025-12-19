# Penalized Complexity (PC) prior for (log) range

Penalized Complexity (PC) prior for (log) range

## Usage

``` r
pclrange(lrange, lam, d = 2, logdens = FALSE)

pcrange(range, lam, d = 2, logdens = FALSE)
```

## Arguments

- lrange:

  numeric with the log of the (practical) range

- lam:

  numeric with the prior parameter

- d:

  integer to specify the domain dimention

- logdens:

  logical indicating if the density is to be returned in the log scale

- range:

  numeric with the of the (practical) range

## Examples

``` r
# P(range < 2.0) = 0.1
 lam <- -log(0.1) * 2.0
 plot(function(x) pcrange(x, lam), 1/100, 10, n = 100)
```
