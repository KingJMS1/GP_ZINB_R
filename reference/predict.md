# Predict at new locations and times

Predicts from fitted Gaussian-process model output at new spatial
locations and times.

## Usage

``` r
predict(X, Ds_new, Dt_new, Vs_new, Vt_new, output)
```

## Arguments

- X:

  Fixed-effect design matrix for the prediction observations.

- Ds_new:

  Distances between new and observed spatial locations.

- Dt_new:

  Distances between new and observed temporal locations.

- Vs_new:

  Spatial random-effect design matrix for prediction observations.

- Vt_new:

  Temporal random-effect design matrix for prediction observations.

- output:

  Fitted model output containing posterior parameter draws.

## Value

A prediction object. This function is currently a placeholder and does
not return calculated predictions.
