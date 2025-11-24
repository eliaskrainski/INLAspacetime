# Helper functions to retrieve the world map, a world polygon, and create grid centers.

Retrieve the map of the countries

## Usage

``` r
worldMap(
  crs = "+proj=moll +units=km",
  scale = "medium",
  returnclass = c("sf", "sv")
)
```

## Arguments

- crs:

  a string with the projection. Default is the Mollweide projection with
  units in kilometers.

- scale:

  The scale of map to return. Please see the help of 'ne_countries'
  function from the 'rnaturalearth' package.

- returnclass:

  A string determining the class of the spatial object to return. Please
  see the help of 'ne_countries' function from the 'rnaturalearth'
  package.

## References

The land and ocean maps are obtained with the 'rnaturalearth' package.
