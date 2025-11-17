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

## Value

Minimum values for l*s, maximum values for l*t
