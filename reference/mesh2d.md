# Illustrative code for building a mesh in 2d domain.

Creates a mesh object. This is just a test code. For efficient, reliable
and general code use the fmesher package.

## Usage

``` r
mesh2d(loc, domain, max.edge, offset, SP = TRUE)
```

## Arguments

- loc:

  a two column matrix with location coordinates.

- domain:

  a two column matrix defining the domain.

- max.edge:

  the maximum edge length.

- offset:

  the length of the outer extension.

- SP:

  logical indicating if the output will include the SpatialPolygons.

## Value

a mesh object.

## Warning

This is just for illustration purposes and one should consider the
efficient function
[`fmesher::fm_mesh_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_mesh_2d.html)
(and other related functions) available a the fmesher package.
