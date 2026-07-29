# Validate Spatial Coordinate Alignment

Confirms that the full coordinate matrix, the baseline-coded spatial
design, and the spatial distance matrix describe the same number of
spatial levels.

## Usage

``` r
validate_spatial_coordinates(coords, Vs, Ds)
```

## Arguments

- coords:

  Spatial coordinate matrix with one row per full spatial level,
  including the baseline level.

- Vs:

  Spatial GP design matrix with the baseline column omitted.

- Ds:

  Spatial distance matrix with one row and column per row of `coords`.

## Value

`TRUE`, invisibly, when the spatial inputs are aligned. Otherwise, the
function stops with an input-specific error message.
