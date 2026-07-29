# Draw posterior predictions without inflation random effects

Draw a zero-inflated negative-binomial response conditional on fixed
effects and count-component spatial and temporal random effects.

## Usage

``` r
estimate_noab(X, alpha, beta, Vs, Vt, c, d, r)
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

- c:

  Spatial random-effect draw for the count component.

- d:

  Temporal random-effect draw for the count component.

- r:

  Negative-binomial dispersion parameter.

## Value

A vector of posterior predictive count draws.
