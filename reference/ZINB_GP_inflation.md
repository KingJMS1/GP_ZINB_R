# ZINB_GP_inflation

Run the ZINB GP model with GP random effects only in the inflation
component

## Usage

``` r
ZINB_GP_inflation(
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
