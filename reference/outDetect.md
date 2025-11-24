# Detect outliers in a time series considering the raw data and a smoothed version of it.

Detect outliers in a time series considering the raw data and a smoothed
version of it.

## Usage

``` r
outDetect(x, weights = NULL, ff = c(7, 7))
```

## Arguments

- x:

  numeric vector

- weights:

  non-increasing numeric vector used as weights for computing a smoothed
  vector as a rooling window average. Default is null and then \\w_j\\
  is proportional to j in the equation in the Details below.

- ff:

  numeric length two vector with the factors used to consider how many
  times the standard deviation one data point is out to be considered as
  an outlier.

## Value

logical vector indicating if the data is an outlier with attributes as
detailed bellow.

- attr(, 'm') is the mean of x.

- attr(, 's') is the standard devation of x.

- attr(, 'ss') is the standard deviation for the smoothed data \\y_t\\
  that is defined as

\\y_t = \sum\_{k=j}^h w_j \* (x\_{t-j}+x\_{t+j})/2\\

Both `s` and `ss` are used to define outliers if

\\\|x_t-m\|/s\>ff_1\\ or \\\|x_t-y_t\|/ss\>ff_2\\

- attr(, 'xs') the smoothed time series \\y_t\\
