# Draw from the full posterior predictive distribution

Draw a zero-inflated negative-binomial response conditional on sampled
fixed effects and spatial and temporal random effects.

## Usage

``` r
estimate(X, alpha, beta, Vs, Vt, a, b, c, d, r)
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

- c:

  Spatial random-effect draw for the count component.

- d:

  Temporal random-effect draw for the count component.

- r:

  Negative-binomial dispersion parameter.

## Value

A vector of posterior predictive count draws.
