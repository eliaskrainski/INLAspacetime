# Internal util functions for polygon properties.

This computes the area of a triangle given its three coordinates.

## Usage

``` r
Heron(x, y)

Area(x, y)

s2trArea(tr, R = 1)

flatArea(tr)

Stiffness(tr)
```

## Arguments

- x, y:

  coordinate vectors.

- tr:

  the triangle coordinates

- R:

  the radius of the spherical domain

## Value

the area of a 2d triangle

the area of a 2d polygon

the area of a triangle in S2

the area of a triangle

the stiffness matrix for a triangle

## Details

Function used internally to compute the area of a triangle.

## Warning

Internal functions, not exported.
