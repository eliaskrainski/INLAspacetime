# Define a regular grid in 'Mollweide' projection, with units in kilometers.

Define a regular grid in 'Mollweide' projection, with units in
kilometers.

## Usage

``` r
world_grid(size = 50, domain)
```

## Arguments

- size:

  the (in kilometers) of the grid cells.

- domain:

  if provided it should be an `sf` or `sfc` object. In this case, the
  grid cells with centers falling inside will be retrieved.

## Value

a 'sf' points object with the centers of a grid set within Earth (and
the supplied domain)
