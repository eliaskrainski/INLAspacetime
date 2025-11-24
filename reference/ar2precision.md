# Precision matrix for a stationary AR2 model.

Creates a precision matrix as a sparse matrix object considering the
specification stated in Details.

## Usage

``` r
ar2precision(n, a)
```

## Arguments

- n:

  integer with the size of the precision matrix.

- a:

  numeric vector with length three with the coefficients.

## Value

sparse matrix.

## Details

Let the second order auto-regression model be defined as

\$\$a_0 x_t + a_1 x\_{t-1} + a_2 x\_{t-2} = w_t, w_t ~ N(0, 1).\$\$

The stationary assumption is to consider the variance of \\x_t\\ to be
the same for all \\t=1,\ldots,n\\. This assumption gives the \\n \times
n\\ symmetric precision matrix \\Q\\ as a sparse matrix with the
following non-zero elements:

\\Q\_{1,1} = Q\_{n,n} = a_0^2\\

\\Q\_{2,2} = Q\_{n-1,n-1} = a_0^2 + a_1^2\\

\\Q\_{1,2} = Q\_{2,1} = Q\_{n-1,n} = Q\_{n,n-1} = a_0 a_1\\

\\Q\_{t,t} = q_0 = a_0^2 + a_1^2 + a_2^2, t = 3, 4, ..., n-2\\

\\Q\_{t,t-1} = Q\_{t-1,t} = q_1 = a_1(a_0 + a_2), t = 3, 4, ..., n-1\\

\\Q\_{t,t-2} = Q\_{t-2,t} = q_2 = a_2 a_0, t = 3, 4, ..., n\\

## See also

[ar2cov](https://eliaskrainski.github.io/INLAspacetime/reference/ar2cov.md)

## Examples

``` r
ar2precision(7, c(1, -1.5, 0.9))
#> 7 x 7 sparse Matrix of class "dgTMatrix"
#>                                             
#> [1,]  1.0 -1.50  0.90  .     .     .     .  
#> [2,] -1.5  3.25 -2.85  0.90  .     .     .  
#> [3,]  0.9 -2.85  4.06 -2.85  0.90  .     .  
#> [4,]  .    0.90 -2.85  4.06 -2.85  0.90  .  
#> [5,]  .    .     0.90 -2.85  4.06 -2.85  0.9
#> [6,]  .    .     .     0.90 -2.85  3.25 -1.5
#> [7,]  .    .     .     .     0.90 -1.50  1.0
```
