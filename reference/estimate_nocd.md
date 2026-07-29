# Draw posterior predictions without count random effects

Draw a zero-inflated negative-binomial response conditional on fixed
effects and zero-inflation-component spatial and temporal random
effects.

## Usage

``` r
estimate_nocd(X, alpha, beta, Vs, Vt, a, b, r)
```

## Arguments

- X:

  Fixed-effect design matrix.

- alpha:

  Fixed-effect coefficients for the zero-inflation component.

- beta:

  Fixed-effect coefficients for the count component.

- Vs:

  Spatial random-effect design matrix.

- Vt:

  Temporal random-effect design matrix.

- a:

  Spatial random-effect draw for the zero-inflation component.

- b:

  Temporal random-effect draw for the zero-inflation component.

- r:

  Negative-binomial dispersion parameter.

## Value

A vector of posterior predictive count draws.
