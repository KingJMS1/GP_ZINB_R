# mvn_sample_svd

Use SVD for precision matrix to sample from multivariate normal
distribution

## Usage

``` r
mvn_sample_svd(P_svd, mu, entropy = NULL, threshold = 1e-12)
```

## Arguments

- P_svd:

  SVD of precision matrix

- mu:

  Mean of MVN to draw from

- entropy:

  Draw from MVN of correct size (can be used to draw all mvns at once
  for efficiency)
