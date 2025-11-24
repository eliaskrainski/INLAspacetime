# Function to define the boundary Earch polygon in longlat projection for a given resolution.

Function to define the boundary Earch polygon in longlat projection for
a given resolution.

## Usage

``` r
Earth_poly(resol = 300, crs = "+proj=moll +units=km")
```

## Arguments

- resol:

  is the number of subdivisions along the latitude coordinates and half
  the number of subdivisions along the longitude coordinates.

- crs:

  a string with the projection. Default is the Mollweide projection with
  units in kilometers.

## Value

a 'st_sfc' object with the Earth polygon.
