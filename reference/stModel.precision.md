# Spacetime precision matrix.

To build the the precision matrix for a spacetime model given the
temporal and the spatial meshes.

## Usage

``` r
stModel.precision(smesh, tmesh, model, theta, verbose = FALSE)
```

## Arguments

- smesh:

  a mesh object over the spatial domain.

- tmesh:

  a mesh object over the time domain.

- model:

  a string identifying the model. So far we have the following models:
  '102', '121', '202' and '220' models.

- theta:

  numeric vector of length three with \\log(gamma_s, gamma_t,
  gamma_e)\\.

- verbose:

  logical to print intermediate objects.

## Value

a (sparse) precision matrix, as in Lindgren et. al. (2024)
