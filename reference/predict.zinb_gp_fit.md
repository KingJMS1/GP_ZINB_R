# Predict at New Locations and Times

Generates posterior predictive draws at new spatial locations, times, or
both. For every retained MCMC iteration, the method draws each active
random effect from its GP predictive distribution conditional on the
stored fitted random effect, adds it to the corresponding fixed-effect
predictor, and draws a zero-inflated negative-binomial response.

## Usage

``` r
# S3 method for class 'zinb_gp_fit'
predict(
  object,
  X,
  Ds_new = NULL,
  Dt_new = NULL,
  Vs_new = NULL,
  Vt_new = NULL,
  kern = NULL,
  ...
)
```

## Arguments

- object:

  A fitted zinb_gp_fit object returned by ZINB_GP().

- X:

  Fixed-effect design matrix for the prediction observations.

- Ds_new:

  Augmented spatial distance matrix produced by
  make_prediction_inputs(), or NULL for a temporal-only fit.

- Dt_new:

  Augmented temporal distance matrix produced by
  make_prediction_inputs(), or NULL for a spatial-only fit.

- Vs_new:

  Spatial design matrix produced by make_prediction_inputs().

- Vt_new:

  Temporal design matrix produced by make_prediction_inputs().

- kern:

  Kernel function used to fit the model. The default is kernel(); supply
  the fitting kernel again if a custom kernel was used.

- ...:

  Unused.

## Value

An object of class zinb_gp_prediction containing posterior draws Y_pred,
eta_at_risk, and eta_count. Active predicted GP effects are returned
under A, B, C, and D, matching the fitted object.

## Examples

``` r
if (FALSE) { # \dontrun{
inputs <- make_prediction_inputs(
  coords = coords,
  time_coords = time_coords,
  coords_new = coords_new,
  time_coords_new = time_coords_new
)
predictions <- predict(
  fit,
  X = X_new,
  Ds_new = inputs$Ds_new,
  Dt_new = inputs$Dt_new,
  Vs_new = inputs$Vs_new,
  Vt_new = inputs$Vt_new
)
} # }
```
