# Define a spacetime model object for the `f()` call.

Define a spacetime model object for the
[`f()`](https://rdrr.io/pkg/INLA/man/f.html) call.

## Usage

``` r
stModel.define(
  smesh,
  tmesh,
  model,
  control.priors,
  constr = FALSE,
  debug = FALSE,
  useINLAprecomp = TRUE,
  libpath = NULL
)
```

## Arguments

- smesh:

  a spatial mesh

- tmesh:

  a temporal mesh

- model:

  a three characters string to specify the smoothness alpha (each one as
  integer) parameters. Currently it considers the `102`, `121`, `202`
  and `220` models.

- control.priors:

  a named list with parameter priors, named as `prs`, `prt` and
  `psigma`, each one as a vector with length two containing (U, a) to
  define the corresponding PC-prior such that, respectively,
  P(range.spatial\<U)=a, P(range.temporal\<U)=a or P(sigma\>U)=a. If a=0
  or a=NA, then U is taken to be the fixed value of the parameter.

- constr:

  logical, default is FALSE, to indicate if the integral of the field
  over the domain is to be constrained to zero.

- debug:

  integer, default is zero, indicating the verbose level. Will be used
  as logical by INLA.

- useINLAprecomp:

  logical, default is TRUE, indicating if it is to be used the shared
  object pre-compiled by INLA. This is not considered if 'libpath' is
  provided.

- libpath:

  string, default is NULL, with the path to the shared object.

## Value

objects to be used in the f() formula term in INLA.

## Details

This function compute the matrices for computing the precision matrix.
These are each one of the Kronecker products in Theorem 4.1 of Lindgren
et. al. (2024) computed with the
[stModel.matrices](https://eliaskrainski.github.io/INLAspacetime/reference/stModel.matrices.md)
and the parameters are as in Eq (19-21). We use the log of these
parameters internally.

## References

Finn Lindgren, Haakon Bakka, David Bolin, Elias Krainski and Håvard Rue
(2024). A diffusion-based spatio-temporal extension of Gaussian Matérn
fields. [SORT vol. 48, no. 1, pp.
3-66](https://raco.cat/index.php/SORT/article/view/428665) \<doi:
10.57645/20.8080.02.13\>
