# make_y_Vs_Vt

Create y, along with spatial and temporal design matrices from an
observation matrix.

## Usage

``` r
make_y_Vs_Vt(obs_matrix)
```

## Arguments

- obs_matrix:

  s by t matrix, where s is the number of locations, t is the number of
  times, each entry of the matrix is a nonnegative integer.

## Value

A list with the following elements:

- y:

  The observation matrix flattened in column-major order.

- Vs:

  Spatial design matrix mapping observations to locations.

- Vt:

  Temporal design matrix mapping observations to times.
