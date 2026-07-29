# Validate a GP Design and Distance Matrix

Checks that a baseline-coded GP design is aligned with the observations
and that its distance matrix contains one additional row and column for
the omitted baseline level.

## Usage

``` r
validate_gp_design(V, D, n_observations, design_name, distance_name)
```

## Arguments

- V:

  GP design matrix with one row per observation and one column per
  nonbaseline GP level.

- D:

  Pairwise distance matrix for all GP levels, including the baseline.

- n_observations:

  Expected number of rows in `V`.

- design_name:

  Character label used to identify `V` in error messages.

- distance_name:

  Character label used to identify `D` in error messages.

## Value

`TRUE`, invisibly, when the design and distance matrix are valid.
Otherwise, the function stops with an input-specific error message.
