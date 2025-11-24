# Functions to help converting from/to user/internal parametrization. The internal parameters are 'gamma_s, 'gamma_t', 'gamma_E' The user parameters are 'r_s', 'r_t', 'sigma'

Functions to help converting from/to user/internal parametrization. The
internal parameters are 'gamma_s, 'gamma_t', 'gamma_E' The user
parameters are 'r_s', 'r_t', 'sigma'

Convert from user parameters to SPDE parameters

Convert from SPDE parameters to user parameters

## Usage

``` r
lgsConstant(lg.s, alpha, smanifold)

params2gammas(
  lparams,
  alpha.t,
  alpha.s,
  alpha.e,
  smanifold = "R2",
  verbose = FALSE
)

gammas2params(lgammas, alpha.t, alpha.s, alpha.e, smanifold = "R2")
```

## Arguments

- lg.s:

  the logarithm of the SPDE parameter `\gamma_s`

- alpha:

  the resulting spatial order.

- smanifold:

  spatial domain manifold, which could be "S1", "S2", "R1", "R2" and
  "R3".

- lparams:

  log(spatial range, temporal range, sigma)

- alpha.t:

  temporal order of the SPDE

- alpha.s:

  spatial order of the spatial differential operator in the
  non-separable part.

- alpha.e:

  spatial order of the spatial differential operator in the separable
  part.

- verbose:

  logical if it is to print internal variables

- lgammas:

  numeric of length 3 with \\log(\gamma_k)\\ model parameters. The
  parameter order is log(gamma.s, gamma.t, gamma.e)

## Value

the part of `sigma` from the spatial constant and `\gamma_s`.

log(gamma.s, gamma.t, gamma.e)

log(spatial range, temporal range, sigma)

## Details

See equation (23) in the paper.

See equations (19), (20) and (21) in the paper.

See equations (19), (20) and (21) in the paper.

## Examples

``` r
params2gammas(log(c(1, 1, 1)), 1, 2, 1, "R2")
#>   gamma.s   gamma.t   gamma.e 
#>  1.039721  1.386294 -3.344954 
gammas2params(log(c(0, 0, 0)), 1, 2, 1, "R2")
#>    lrs    lrt lsigma 
#>    Inf    NaN    Inf 
```
