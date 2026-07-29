# Draw a Zero-Inflated Negative-Binomial Response

Generates one posterior predictive count for each pair of at-risk and
count linear predictors using the parameterization fitted by
[`ZINB_GP()`](https://kingjms1.github.io/GP_ZINB_R/reference/ZINB_GP.md).

## Usage

``` r
draw_zinb_response(eta_at_risk, eta_count, r)
```

## Arguments

- eta_at_risk:

  Numeric vector of linear predictors for the Bernoulli at-risk
  component.

- eta_count:

  Numeric vector of linear predictors for the negative-binomial count
  component.

- r:

  Positive negative-binomial dispersion parameter.

## Value

A non-negative integer vector with one predictive count per linear
predictor pair.
