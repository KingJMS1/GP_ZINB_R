# ZINB (NN, in future)GP Bayesian Model

This package implements the model first described in “A Framework of
Zero-Inflated Bayesian Negative Binomial Regression Models For
Spatiotemporal Data” by Qing He and Hsin-Hsiung Huang (2023).
<https://doi.org/10.1016/j.jspi.2023.106098>, later amended by
“Imaging-derived spatiotemporal modeling of landslide counts with noisy
Gaussian process random effects.”

This package is a work in progress, feel free to create an issue if you
have suggestions or notice any problems.

## Installation Instructions

Install the devtools R package, then run the following command

``` r
install.packages(c("BayesLogit", "LaplacesDemon", "MASS", "Matrix", "msm", "mvtnorm"))
devtools::install_github("KingJMS1/NNGP_ZINB_R")
```

## Example Use

Detailed examples with full code can be found in the vignettes folder,
source and data for which can be found at the experiments repository:
<https://github.com/KingJMS1/ZINB_R_Experiments>.

    X          Other Predictor variables
    y          Zero inflated count response
    coords     Spatial coordinates for GP
    Vs         Spatially varying predictor variables
               (e.g. one-hot indication of which location this is for varying intercept),
               wrapped in sparseMatrix from Matrix R package.
               Will be multiplied by the spatial random effects for prediction.
    Vt         Temporal varying predictor variables, wrapped in sparseMatrix from Matrix R package.
               Will be multiplied by the temporal random effects for prediction.
    Ds         Spatial distance matrix, diagonal should be 0,
               off diagonal is distance between elements i and j in space, inputs to the spatial GP kernel
    Dt         Temporal distance matirx, diagonal should be 0,
               off diagonal is distance between elements i and j in time, inputs to the temporal GP kernel
    nsim       How long to run MCMC in total, must be greater than burn.
    burn       How long to run MCMC before saving samples.
    thin       Thinning ratio for chain
    ltPrior    Prior for temporal length scale, formatted as: list(max=50, mh_sd=3, a=1, b=0.001). a/b for gamma prior, max sets how high it is allowed to go, mh_sd controls the proposal variance for mh.
    lsPrior    Prior for spatial length scale, formatted the same as ltPrior.
    sigmaPrior Prior for sigmas in model, formatted as: list(a=0.01, b=0.1), inverse gamma prior
    noisePrior Prior for signal to noise ratio, formatter as: list(a=1.5, b=1.5, mh_sd=0.2), mh_sd controls the proposal variance for mh.
    mh_sd_r    Control mh proposal variance for r to deal with convergence issues
    kern       Covariance kernel, takes a squared distance matrix, returns kernel covariance.

Given all of the above, predictions can then be found via:

``` r
library(ZINB.GP)
output <- ZINB_GP(X, y, coords, Vs, Vt, Ds, Dt, M = M, nsim, burn, save_ypred = TRUE)
predictions <- output$Y_pred
```

## API Reference

API reference:
<https://kingjms1.github.io/NNGP_ZINB_R/reference/index.html>
