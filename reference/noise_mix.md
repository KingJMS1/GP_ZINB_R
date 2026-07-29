# Mix a Matrix with Identity Noise

Creates `noise_ratio * A + (1 - noise_ratio) * I`.

## Usage

``` r
noise_mix(A, noise_ratio)
```

## Arguments

- A:

  Numeric square matrix.

- noise_ratio:

  Numeric mixing ratio between zero and one.

## Value

A numeric square matrix with the same dimensions as `A`.

## Examples

``` r
correlation <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
noise_mix(correlation, noise_ratio = 0.8)
#>      [,1] [,2]
#> [1,]  1.0  0.4
#> [2,]  0.4  1.0
```
