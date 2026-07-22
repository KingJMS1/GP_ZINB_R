# ZINB_GP_spatial2

Fit a zero-inflated negative-binomial model with spatial
Gaussian-process random effects in only the count component.

## Usage

``` r
ZINB_GP_spatial2(
  X,
  y,
  Vs,
  Ds,
  nsim,
  burn,
  thin = 1,
  save_ypred = FALSE,
  print_iter = 100,
  print_progress = FALSE,
  lsPrior = NULL,
  sigmaPrior = NULL,
  noisePrior = NULL,
  mh_sd_r = NULL,
  kern = NULL
)
```

## Arguments

- X:

  Fixed-effect design matrix with N rows.

- y:

  Non-negative integer count response of length N.

- Vs:

  Sparse or dense spatial random-effect design matrix. It should have N
  rows and one column per spatial location.

- Ds:

  Spatial distance matrix, with one row and column per spatial random
  effect. Diagonal entries must be zero.

- nsim:

  Total number of MCMC iterations.

- burn:

  Number of burn-in iterations.

- thin:

  Store every thin-th iteration after burn-in.

- save_ypred:

  Whether to save posterior fitted means and at-risk draws.

- print_iter:

  Print progress every print_iter iterations.

- print_progress:

  Whether to print MCMC progress.

- lsPrior:

  Prior and proposal controls for spatial GP length scales.

- sigmaPrior:

  Inverse-gamma prior parameters for GP variances.

- noisePrior:

  Beta prior and MH controls for GP noise ratios.

- mh_sd_r:

  Proposal standard deviation for NB dispersion r.

- kern:

  Kernel function accepting a squared-distance matrix and a length
  scale.

## Value

A list containing posterior MCMC draws.
