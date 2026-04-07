# Build a `cgeneric` for Nearest Neighbour Gaussian Process.

Build data needed to implement the model for the C interface,
`cgeneric`, in the 'INLA' package, that can be used as a linear
predictor model component with an 'f' term.

## Usage

``` r
cgeneric_nngp(dists, nnl, control.priors, correl = "Matern", nu = 1, ...)
```

## Arguments

- dists:

  a `dist` object.

- nnl:

  a list with the neighbors indexes.

- control.priors:

  a named list with parameter priors, named as `prs`, `prt` and
  `psigma`, each one as a vector with length two containing (U, a) to
  define the corresponding PC-prior such that, respectively,
  P(range.spatial\<U)=a, P(rho\>U)=a or P(sigma\>U)=a. If a=0 or a=NA,
  then U is taken to be the fixed value of the parameter.

- correl:

  character to specify the correlation function. Available options are
  "Matern" (default) and "powerExponential".

- nu:

  numeric to specify the smoothness parameter.

- ...:

  arguments (debug,useINLAprecomp,shlib) passed on to
  [`INLAtools::cgeneric()`](https://rdrr.io/pkg/INLAtools/man/cgeneric-class.html).

## Value

a `cgeneric` object, see
[`INLAtools::cgeneric()`](https://rdrr.io/pkg/INLAtools/man/cgeneric-class.html).

## Details

The precision matrix is defined as \$\$Q = LL',L = (D-A)D^{-0.5}\$\$.
where A is a lower triangle matrix and D a diagonal one..

## See also

[`INLAtools::prior.cgeneric()`](https://rdrr.io/pkg/INLAtools/man/cgeneric_get.html)
