# Compute LBM covariance matrix from distance matrix

sigma_lbm = 0.5 \* (J \* dist,1' + dist1, \* J' - dist) Uses the first
location as the origin.

## Usage

``` r
sigma_lbm(distmat)
```

## Arguments

- distmat:

  n x n distance matrix

## Value

n x n LBM covariance matrix
