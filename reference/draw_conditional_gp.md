# Draw a GP at New Levels Conditional on a Stored Draw

Draw a GP at New Levels Conditional on a Stored Draw

## Usage

``` r
draw_conditional_gp(effect, distance, length_scale, sigma, noise_ratio, kern)
```

## Arguments

- effect:

  Numeric GP draw at the observed nonbaseline levels.

- distance:

  Augmented distance matrix whose observed levels precede its new
  levels.

- length_scale:

  GP length scale.

- sigma:

  GP marginal standard deviation.

- noise_ratio:

  Fraction of variance assigned to the structured kernel.

- kern:

  Kernel function accepting squared distances and a length scale.

## Value

A numeric draw at the new levels.
