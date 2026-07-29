# ZINB_GP_count

Run the ZINB GP model with GP random effects only in the count part.

## Usage

``` r
ZINB_GP_count(
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
