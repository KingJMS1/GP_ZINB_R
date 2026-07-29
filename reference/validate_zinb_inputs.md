# Validate ZINB Model Inputs

Checks the fixed-effect design, count response, and MCMC controls shared
by the branch-specific ZINB-GP samplers.

## Usage

``` r
validate_zinb_inputs(X, y, nsim, burn, thin)
```

## Arguments

- X:

  Fixed-effect design matrix with one row per observation.

- y:

  Non-negative integer count response with one value per row of `X`.

- nsim:

  Positive integer giving the total number of MCMC iterations.

- burn:

  Non-negative integer giving the number of burn-in iterations.

- thin:

  Positive integer giving the post-burn-in thinning interval.

## Value

`TRUE`, invisibly, when all inputs are valid. Otherwise, the function
stops with an input-specific error message.
