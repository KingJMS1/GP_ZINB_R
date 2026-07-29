# make_y_Vs_Vt

Create y, along with spatial and temporal design matrices from an
observation matrix.

## Usage

``` r
make_y_Vs_Vt(obs_matrix)
```

## Arguments

- obs_matrix:

  s by t matrix, where s is the number of locations, t is the number of
  times, each entry of the matrix is a nonnegative integer.

## Value

A list with the following elements:

- y:

  The observation matrix flattened in column-major order.

- Vs:

  Spatial design matrix mapping observations to locations.

- Vt:

  Temporal design matrix mapping observations to times.

## Examples

``` r
counts <- matrix(c(0, 1, 2, 0, 1, 3), nrow = 2)
design <- make_y_Vs_Vt(counts)
design$y
#> [1] 0 1 2 0 1 3
design$Vs
#> [1] 0 1 0 1 0 1
design$Vt
#>      [,1] [,2]
#> [1,]    0    0
#> [2,]    0    0
#> [3,]    1    0
#> [4,]    1    0
#> [5,]    0    1
#> [6,]    0    1
```
