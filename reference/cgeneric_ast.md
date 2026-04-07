# Build an `cgeneric` conditional autoregression spacetime model.

Build data needed to implement the model for the C interface,
`cgeneric`, in the 'INLA' package, that can be used as a linear
predictor model component with an 'f' term.

## Usage

``` r
cgeneric_ast(G, nt, control.priors, ...)
```

## Arguments

- nt:

  integer to define the temporal dimention.

- control.priors:

  a named list with parameter priors, named as `prs`, `prt` and
  `psigma`, each one as a vector with length two containing (U, a) to
  define the corresponding PC-prior such that, respectively,
  P(range.spatial\<U)=a, P(rho\>U)=a or P(sigma\>U)=a. If a=0 or a=NA,
  then U is taken to be the fixed value of the parameter.

- ...:

  arguments (debug,useINLAprecomp,shlib) passed on to
  [`INLAtools::cgeneric()`](https://rdrr.io/pkg/INLAtools/man/cgeneric-class.html).

- R:

  a graph representing the spatial neighborhood.

## Value

a `cgeneric` object, see
[`INLAtools::cgeneric()`](https://rdrr.io/pkg/INLAtools/man/cgeneric-class.html).

## Details

The precision matrix is defined as \$\$Q = \tau ((D(\rho) + H(\rho))
\otimes (\kappa^2 I + R)\$\$ where D and H are temporal matrices, and R
is the structure matrix supplied. \\\tau\\ is the (local) precision
parameter, \\\kappa^2\\ is a scaling parameter and \\\rho\\ is the
temporal autoregressive parameter.

## See also

[`INLAtools::prior.cgeneric()`](https://rdrr.io/pkg/INLAtools/man/cgeneric_get.html)
