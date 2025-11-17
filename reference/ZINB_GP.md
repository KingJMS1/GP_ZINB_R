# ZINB_GP

Run the ZINB NNGP model described in
https://doi.org/10.1016/j.jspi.2023.106098.

## Usage

``` r
ZINB_GP(
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

  Other Predictor variables

- y:

  Zero inflated count response

- coords:

  Spatial coordinates for NNGP

- Vs:

  Spatially varying predictor variables (e.g. one-hot indication of
  which location this is for varying intercept), wrapped in sparseMatrix
  from Matrix R package. Will be multiplied by the spatial random
  effects for prediction.

- Vt:

  Temporal varying predictor variables, wrapped in sparseMatrix from
  Matrix R package. Will be multiplied by the temporal random effects
  for prediction.

- Ds:

  Spatial distance matrix, diagonal should be 0, off diagonal is
  distance between elements i and j in space, inputs to the spatial NNGP
  kernel

- Dt:

  Temporal distance matirx, diagonal should be 0, off diagonal is
  distance between elements i and j in time, inputs to the temporal GP
  kernel

- nsim:

  How long to run MCMC in total, must be greater than burn.

- burn:

  How long to run MCMC before saving samples.

- thin:

  How often to save MCMC samples, default is 1, saves every iteration.

- save_ypred:

  Whether or not to output the predicted values at every iteration

- print_iter:

  How often to print the iteration number of the MCMC chain.

- print_progress:

  Whether or not to print the iteration number of the MCMC chain.

- M:

  How many neighbors to allow in the spatial NNGP algorithm, defaults to
  10.

## Value

A List of the following sampled values:

- **Alpha:** Model coefficients for logit model

- **Beta:** Model coefficients for NB model

- **A:** Portion of spatial random effect in the logit model explained
  by kernel

- **B:** Portion of temporal random effect in the logit model explained
  by kernel

- **C:** Portion of spatial random effect in the NB model explained by
  kernel

- **D:** Portion of temporal random effect in the NB model explained by
  kernel

- **L1t:** Length scale for temporal kernel in logit model, i.e.
  \\e^{-\frac{d^{2}}{2 l\_{1t}^{2}}}\\

- **Sigma1t:** Kernel scale parameter for above kernel, i.e.
  \\\sigma\_{1t}^{2}e^{.}\\

- **L2t:** Length scale for temporal kernel in NB model, i.e.
  \\e^{-\frac{d^{2}}{2 l\_{1t}^{2}}}\\

- **Sigma2t:** Kernel scale parameter for above kernel, i.e.
  \\\sigma\_{2t}^{2}e^{.}\\

- **Phi_bin:** Length scale for spatial kernel in logit model, i.e.
  \\e^{-\Phi\_{bin}d^{2}}\\

- **Sigma1s:** Square root of multiplier for spatial kernel in logit
  model

- **Phi_nb:** Length scale for spatial kernel in NB model, i.e.
  \\e^{-\Phi\_{nb}d^{2}}\\

- **Sigma2s:** Square root of multiplier for spatial kernel in NB model

- **R:** Dispersion parameter for Negative Binomial distribution.

- **at_risk:** At risk indicator for each observation

- **Y_pred:** Predictions, sampled from the posterior distribution at
  each iteration, NULL if save_ypred is false
