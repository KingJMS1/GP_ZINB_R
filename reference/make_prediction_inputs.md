# Build Inputs for Prediction at New GP Coordinates

Creates the augmented distance matrices needed for Gaussian-process
conditioning and indicator matrices mapping prediction observations to
unique new spatial and temporal levels. The first observed level in each
dimension is the baseline omitted by
[`ZINB_GP()`](https://kingjms1.github.io/GP_ZINB_R/reference/ZINB_GP.md)
and is therefore excluded from the augmented distance matrix.

## Usage

``` r
make_prediction_inputs(
  coords = NULL,
  time_coords = NULL,
  coords_new = NULL,
  time_coords_new = NULL
)
```

## Arguments

- coords:

  Observed spatial coordinate matrix, including the baseline as its
  first row. Use `NULL` for a temporal-only model.

- time_coords:

  Observed temporal coordinate matrix, including the baseline as its
  first row. A numeric vector is treated as a one-column matrix. Use
  `NULL` for a spatial-only model.

- coords_new:

  Spatial coordinates for the prediction observations, with one row per
  row of the new fixed-effect design matrix. Repeated rows share one
  predicted spatial random effect.

- time_coords_new:

  Temporal coordinates for the prediction observations, with one row per
  row of the new fixed-effect design matrix. A numeric vector is treated
  as a one-column matrix. Repeated rows share one predicted temporal
  random effect.

## Value

A list containing `Ds_new`, `Dt_new`, `Vs_new`, and `Vt_new`. Each
non-`NULL` distance matrix covers the observed nonbaseline levels
followed by the unique new levels. Each design matrix has one row per
prediction observation and one column per unique new level.

## Examples

``` r
observed_space <- rbind(c(0, 0), c(1000, 0), c(0, 1000))
observed_time <- c(0, 100)
new_space <- rbind(c(500, 500), c(500, 500), c(750, 250))
new_time <- c(150, 200, 150)
prediction_inputs <- make_prediction_inputs(
  coords = observed_space,
  time_coords = observed_time,
  coords_new = new_space,
  time_coords_new = new_time
)
prediction_inputs$Vs_new
#>      [,1] [,2]
#> [1,]    1    0
#> [2,]    1    0
#> [3,]    0    1
prediction_inputs$Vt_new
#>      [,1] [,2]
#> [1,]    1    0
#> [2,]    0    1
#> [3,]    1    0
```
