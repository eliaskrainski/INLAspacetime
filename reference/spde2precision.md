# Illustrative code to build the precision matrix for SPDE kind models.

Creates a precision matrix as a sparse matrix object. For general code
look at the functions in the INLA package.

## Usage

``` r
spde2precision(kappa, fem, alpha)
```

## Arguments

- kappa:

  the scale parameter.

- fem:

  a list containing the Finite Element matrices.

- alpha:

  the smoothness parameter.

## Value

the precision matrix as a sparse matrix object.

## Warning

This is just for illustration purpose and one should consider the
efficient function available a the INLA package.
