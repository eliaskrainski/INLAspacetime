# Illustrative code to build the projector matrix for SPDE models.

Creates a projector matrix object.

## Usage

``` r
mesh2projector(
  mesh,
  loc = NULL,
  lattice = NULL,
  xlim = NULL,
  ylim = NULL,
  dims = c(100, 100)
)
```

## Arguments

- mesh:

  a 2d mesh object.

- loc:

  a two columns matrix with the locations to project for.

- lattice:

  Unused; feature not supported by this illustration.

- xlim, ylim:

  vector with the boundary limits.

- dims:

  the number of subdivisions over each boundary limits.

## Value

the projector matrix as a list with sparse matrix object at `x$proj$A`..

## Warning

This is just for illustration purpose and one should consider the
efficient functions available in the fmesher package, e.g.
[`fmesher::fm_evaluator()`](https://inlabru-org.github.io/fmesher/reference/fm_evaluate.html)
and
[`fmesher::fm_basis()`](https://inlabru-org.github.io/fmesher/reference/fm_basis.html).
