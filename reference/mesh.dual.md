# Extracts the dual of a mesh object.

Extracts the dual of a mesh object.

## Usage

``` r
mesh.dual(
  mesh,
  returnclass = c("list", "sf", "sv", "SP"),
  mc.cores = getOption("mc.cores", 2L)
)
```

## Arguments

- mesh:

  a 2d mesh object.

- returnclass:

  if 'list' return a list of polygon coordinates, if "sf" return a 'sf'
  sfc_multipolygon object, if "sv" return a 'terra', SpatVector object,
  if "SP" return a 'sp' SpatialPolygons object.

- mc.cores:

  number of threads to be used.

## Value

one of the three in 'returnclass'
