# ZINB_GP

Fits the zero-inflated negative-binomial Gaussian-process model
described in <https://doi.org/10.1016/j.jspi.2023.106098>. The supplied
spatial and temporal design and distance matrices determine which
Gaussian processes are included.

## Usage

``` r
ZINB_GP(
  X,
  y,
  coords,
  nsim = 5000,
  burn = 1000,
  use_count_gp = TRUE,
  use_inflation_gp = FALSE,
  thin = 1,
  kern = NULL,
  save_ypred = FALSE,
  print_iter = 100,
  print_progress = FALSE,
  Vs = NULL,
  Vt = NULL,
  Ds = NULL,
  Dt = NULL,
  ltPrior = NULL,
  lsPrior = NULL,
  sigmaPrior = NULL,
  noisePrior = NULL,
  mh_sd_r = NULL
)
```

## Arguments

- X:

  Fixed-effect design matrix with one row per observation.

- y:

  Non-negative integer count response.

- coords:

  Spatial coordinates used by the legacy spNNGP implementation.

- nsim:

  Total number of MCMC iterations; must exceed `burn`.

- burn:

  Number of burn-in iterations.

- use_count_gp:

  Whether to include GP random effects in the count component.

- use_inflation_gp:

  Whether to include GP random effects in the zero-inflation component.

- thin:

  Store every `thin`-th iteration after burn-in.

- kern:

  Kernel function accepting a distance matrix and length scale.

- save_ypred:

  Whether to save posterior predictive draws.

- print_iter:

  Print progress every `print_iter` iterations.

- print_progress:

  Whether to print MCMC progress.

- Vs:

  Spatial random-effect design matrix with one row per observation.

- Vt:

  Temporal random-effect design matrix with one row per observation.

- Ds:

  Spatial distance matrix; its diagonal must be zero.

- Dt:

  Temporal distance matrix; its diagonal must be zero.

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

## Value

A list containing posterior MCMC draws:

- Alpha:

  Fixed-effect coefficients for the zero-inflation component.

- Beta:

  Fixed-effect coefficients for the count component.

- A, B:

  Spatial and temporal random effects for the zero-inflation component.

- C, D:

  Spatial and temporal random effects for the count component.

- L1t, L2t:

  Temporal GP length scales for the zero-inflation and count components.

- Sigma1t, Sigma2t:

  Temporal GP scale parameters.

- Phi_bin, Phi_nb:

  Spatial GP length-scale parameters.

- Sigma1s, Sigma2s:

  Spatial GP scale parameters.

- R:

  Negative-binomial dispersion parameter.

- at_risk:

  At-risk indicator draws for each observation.

- Y_pred:

  Posterior predictive draws, or `NULL` when `save_ypred` is `FALSE`.
