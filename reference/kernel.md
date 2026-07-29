# Squared Exponential Kernel

Creates the kernel matrix `exp(-dist / ls^2)`.

## Usage

``` r
kernel(dist, ls)
```

## Arguments

- dist:

  Numeric distance matrix.

- ls:

  Positive length scale.

## Value

A numeric kernel matrix with the same dimensions as `dist`.

## Examples

``` r
distances <- as.matrix(stats::dist(c(0, 1, 3)))
kernel(distances, ls = 2)
#>           1         2         3
#> 1 1.0000000 0.7788008 0.4723666
#> 2 0.7788008 1.0000000 0.6065307
#> 3 0.4723666 0.6065307 1.0000000
```
