# gp_param_bounds

Finds reasonable upper/lower bounds on gp parameters ensuring matrices
remain pd invertible

## Usage

``` r
gp_param_bounds(Ds, Dt, kernel, tolerance = 1e-10)
```

## Arguments

- Ds:

  Spatial distance matrix

- Dt:

  Temporal distance matrix

- kernel:

  Kernel function accepting a distance matrix and a length scale.

- tolerance:

  Numerical tolerance used to assess matrix invertibility.

## Value

A list containing the upper bounds `ltmax` and `lsmax` for the temporal
and spatial length scales, respectively.
