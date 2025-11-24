# Mapper object for automatic inlabru interface

Return an `inlabru` `bru_mapper` object that can be used for computing
model matrices for the space-time model components. The
`bru_get_mapper()` function is called by the `inlabru` methods to
automatically obtain the needed mapper object (from `inlabru`
`2.7.0.9001`; before that, use `mapper = bru_get_mapper(model)`
explicitly).

## Usage

``` r
# S3 method for class 'stModel_cgeneric'
bru_get_mapper(model, ...)
```

## Arguments

- model:

  The model object (of class `stModel_cgeneric`, from `stModel.define`
  or `barrierModel_cgeneric`, from `barrierModel.define`)

- ...:

  Unused.

## Value

A `bru_mapper` object of class `bru_mapper_multi` with sub-mappers
`space` and `time` based on the model `smesh` and `tmesh` or `mesh`
objects.

## See also

[`inlabru::bru_get_mapper()`](https://inlabru-org.github.io/inlabru/reference/bru_get_mapper.html)
