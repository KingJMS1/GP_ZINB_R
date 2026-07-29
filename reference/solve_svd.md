# solve_svd

Use SVD to solve a linear system Ax=b

## Usage

``` r
solve_svd(A_svd, b, threshold = 1e-12)
```

## Arguments

- A_svd:

  SVD of A

- b:

  Right-hand-side vector or matrix.

- threshold:

  Singular values at or below this threshold are discarded.

## Value

The SVD-based solution to `A %*% x = b`.
