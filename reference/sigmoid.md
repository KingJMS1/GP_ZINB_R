# Numerically stable sigmoid function

Computes the inverse-logit transformation while clipping its output away
from zero and one.

## Usage

``` r
sigmoid(eta)
```

## Arguments

- eta:

  Numeric vector or matrix of linear predictors.

## Value

`eta` transformed to probabilities in `[1e-6, 1 - 1e-6]`.
