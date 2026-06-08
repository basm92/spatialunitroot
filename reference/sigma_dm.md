# Compute demeaned spatial covariance with exponential decay

sigma = exp(-c \* distmat), then double-center.

## Usage

``` r
sigma_dm(distmat, c)
```

## Arguments

- distmat:

  n x n distance matrix

- c:

  Decay parameter

## Value

Double-centered exponential covariance matrix
