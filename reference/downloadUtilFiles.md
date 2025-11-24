# Download files from the NOAA's GHCN daily data

Download files from the NOAA's GHCN daily data

## Usage

``` r
downloadUtilFiles(data.dir, year = 2022, force = FALSE)
```

## Arguments

- data.dir:

  the folder to store the files.

- year:

  the year of the daily weather data.

- force:

  logical indicating if it is to force the download. If FALSE each file
  will be downloaded if it does not exists locally yet.

## Value

a named character vector with the local file names: daily.data,
stations.all, elevation.

## See also

[`ghcndSelect()`](https://eliaskrainski.github.io/INLAspacetime/reference/ghcndSelect.md)
