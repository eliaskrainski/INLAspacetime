# To check unusual low/high variance segments

To check unusual low/high variance segments

## Usage

``` r
stdSubs(x, nsub = 12, fs = 15)
```

## Arguments

- x:

  numeric vector

- nsub:

  number for the segments length

- fs:

  numeric to use for detecting too hight or too low local standard
  deviations.

## Value

logical indicating if any of the `st` are `fs` times lower/higher the
average of `st`, where is returned as an attribute:

- attr(, 'st') numeric vector with the standard deviation at each
  segment of the data.
