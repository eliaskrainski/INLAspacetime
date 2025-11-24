# Define the spacetime model matrices.

This function computes all the matrices needed to build the precision
matrix for spatio-temporal model, as in Lindgren et. al. (2024)

## Usage

``` r
stModel.matrices(smesh, tmesh, model, constr = FALSE)
```

## Arguments

- smesh:

  a mesh object over the spatial domain.

- tmesh:

  a mesh object over the time domain.

- model:

  a string identifying the model. So far we have the following models:
  '102', '121', '202' and '220' models.

- constr:

  logical to indicate if the integral of the field over the domain is to
  be constrained to zero. Default value is FALSE.

## Value

a list containing needed objects for model definition.

1.  'manifold' to spedify the a string with the model identification

2.  a length three vector with the constants `c1`, `c2` and `c3`

3.  the vector `d`

4.  the matrix `T`

5.  the model matrices `M_1`, ..., `M_m`

## Details

See the paper for details.
