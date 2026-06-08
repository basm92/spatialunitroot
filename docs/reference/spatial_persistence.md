# Compute spatial half-life confidence interval

Implements the Bayesian posterior-based confidence set for the spatial
half-life parameter. Matches Stata spatial_persistence().

## Usage

``` r
spatial_persistence(Z, distmat, emat, level)
```

## Arguments

- Z:

  n x 1 vector of observations

- distmat:

  n x n normalized distance matrix

- emat:

  q x nrep matrix of standard normal draws

- level:

  Confidence level (e.g., 0.95)

## Value

List with ci_lower, ci_upper
