# Define a spacetime model object for the `f()` call.

Define a spacetime model object for the
[`f()`](https://rdrr.io/pkg/INLA/man/f.html) call.

## Usage

``` r
barrierModel.define(
  mesh,
  barrier.triangles,
  prior.range,
  prior.sigma,
  range.fraction = 0.1,
  constr = FALSE,
  debug = FALSE,
  useINLAprecomp = TRUE,
  libpath = NULL
)
```

## Arguments

- mesh:

  a spatial mesh

- barrier.triangles:

  a integer vector to specify which triangles centers are in the barrier
  domain, or a list with integer vector if more than one.

- prior.range:

  numeric vector containing U and a to define the probability statements
  P(range \< U) = a used to setup the PC-prior for range. If a = 0 or a
  = NA, then U is taken to be the fixed value for the range.

- prior.sigma:

  numeric vector containing U and a to define the probability statements
  P(range \> U) = a used to setup the PC-prior for sigma. If a = 0 or a
  = NA, then U is taken to be the fixed value for sigma.

- range.fraction:

  numeric to specify the fraction of the range for the barrier domain.
  Default value is 0.1. This has to be specified with care in order to
  have it small enough to make it act as barrier but not too small in
  order to prevent numerical issues.

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

See the paper.
