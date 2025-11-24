# To visualize time series over space.

To visualize time series over space.

## Usage

``` r
stlines(
  stdata,
  spatial,
  group = NULL,
  nmax.group = NULL,
  xscale = 1,
  yscale = 1,
  colour = NULL,
  ...
)

stpoints(
  stdata,
  spatial,
  group = NULL,
  nmax.group = NULL,
  xscale = 1,
  yscale = 1,
  colour = NULL,
  ...
)
```

## Arguments

- stdata:

  matrix with the data, each column is a location.

- spatial:

  an object with one of class defined in the sp package.

- group:

  an integer vector indicating to which spatial unit each time series
  belongs to. Default is NULL and them it is assumed that each time
  series belongs o each spatial unit.

- nmax.group:

  an integer indicating the maximum number of time series to be plotted
  over each spatial unit. Default is NULL, so all will be drawn.

- xscale:

  numeric to define a scaling factor in the horizontal direction.

- yscale:

  numeric to define a scaling factor in the vertical direction.

- colour:

  color (may be a vector, one for each time series). Default is NULL and
  it will generate colors considering the average of each time series.
  These automatic colors are defined using the
  [`rgb()`](https://rdrr.io/r/grDevices/rgb.html) function with
  `alpha=0.5`. It considers the relative rank of each time series mean,
  `r`. `r` is then used for red, `1-r` is used for blue and a triangular
  function, `1-2*|1-r/2|`, is considered for green. That is, time series
  with mean among the lowest time series averages are shown in blue and
  those among the highest temperatures are shown in red. The transition
  from blue to red goes so that the intermediate ones are shown in light
  green.

- ...:

  further arguments to be passed for the lines function.

## Value

add lines to an existing plot

## Details

Scaling the times series is needed before drawing it over the map. The
area of the bounding box for the spatial object divided by the number of
locations is the standard scaling factor. This is further multiplied by
the user given `xcale` and `yscale`.

## Functions

- `stlines()`: each time series over the map centered at the location.

- `stpoints()`: each time series over the map centered at the location.

## Warning

if there are too many geographical locations, it will not look good
