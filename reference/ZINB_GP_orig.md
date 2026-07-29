# ZINB_GP_orig

Fits the spatial-temporal zero-inflated negative-binomial
Gaussian-process model described in
[doi:10.1016/j.jspi.2023.106098](https://doi.org/10.1016/j.jspi.2023.106098)
.

## Usage

``` r
ZINB_GP_orig(
  X,
  y,
  coords,
  Vs,
  Vt,
  Ds,
  Dt,
  nsim,
  burn,
  thin = 1,
  save_ypred = FALSE,
  print_iter = 100,
  print_progress = FALSE,
  ltPrior = NULL,
  lsPrior = NULL,
  sigmaPrior = NULL,
  noisePrior = NULL,
  mh_sd_r = NULL,
  kern = NULL
)
```

## Arguments

- X:

  Fixed-effect design matrix with one row per observation.

- y:

  Non-negative integer count response.

- coords:

  Spatial coordinate matrix with one row per full spatial level,
  including the baseline level omitted from `Vs`.

- Vs:

  Spatial random-effect design matrix with one row per observation.

- Vt:

  Temporal random-effect design matrix with one row per observation.

- Ds:

  Spatial distance matrix; its diagonal must be zero.

- Dt:

  Temporal distance matrix; its diagonal must be zero.

- nsim:

  Total number of MCMC iterations; must exceed `burn`.

- burn:

  Number of burn-in iterations.

- thin:

  Store every `thin`-th iteration after burn-in.

- save_ypred:

  Whether to save posterior predictive draws.

- print_iter:

  Print progress every `print_iter` iterations.

- print_progress:

  Whether to print MCMC progress.

- ltPrior:

  List with `max`, `mh_sd`, `a`, and `b` for temporal length-scale prior
  and proposal controls.

- lsPrior:

  List with `max`, `mh_sd`, `a`, and `b` for spatial length-scale prior
  and proposal controls.

- sigmaPrior:

  List with `a` and `b` inverse-gamma prior parameters for GP scales.

- noisePrior:

  List with `a`, `b`, and `mh_sd` for the GP noise-ratio prior and
  proposal.

- mh_sd_r:

  Proposal standard deviation for the negative-binomial dispersion
  parameter.

- kern:

  Kernel function accepting a distance matrix and length scale.

## Value

A list containing posterior MCMC draws, including fixed-effect
coefficients, spatial and temporal random effects, their GP parameters,
the negative-binomial dispersion parameter, and optional predictive
draws.
