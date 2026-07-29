# Conditional At-Risk Probability for an Observed Zero

Computes the conditional probability that an observation belongs to the
negative-binomial, or at-risk, component given that its observed count
is zero. The calculation is performed on the log-odds scale for
numerical stability.

## Usage

``` r
zero_at_risk_probability(pi, q, r)
```

## Arguments

- pi:

  Prior at-risk probabilities.

- q:

  Negative-binomial zero-probability bases, equal to
  `1 - sigmoid(eta_count)`.

- r:

  Positive negative-binomial dispersion parameter.

## Value

A numeric vector or matrix with the same shape as `pi` and `q`,
containing probabilities between zero and one.
