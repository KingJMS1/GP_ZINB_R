# update_ls_sigma_noise

Update kernel parameters for a GP

## Usage

``` r
update_ls_sigma_noise(
  ls,
  sigma,
  noise_ratio,
  gpdraw,
  K,
  D,
  lsPrior,
  sigmaPrior,
  noisePrior,
  kern
)
```

## Arguments

- ls:

  Current length scale

- sigma:

  Current sigma

- noise_ratio:

  Current noise ratio

- gpdraw:

  Last draw from the gp with these parameters

- K:

  Current kernel matrix

- D:

  Distance matrix

- lsPrior:

  prior information for length scale, needs mh_sd, max, a, b

- sigmaPrior:

  prior information for sigma, needs a, b

- noisePrior:

  prior information for noise_ratio, needs mh_sd, a, b

- kern:

  Kernel function accepting a distance matrix and a length scale.

## Value

A list with the following elements:

- ls:

  Updated length scale.

- sigma:

  Updated GP scale.

- noise_ratio:

  Updated kernel-to-noise mixing ratio.

- K:

  Updated kernel matrix.

- K_inv:

  Inverse of the updated kernel matrix.
