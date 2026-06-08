# Compute residual-projected spatial covariance

sigma = exp(-c \* distmat), then sigma_dm = M \* sigma \* M'.

## Usage

``` r
sigma_residual(distmat, c, M)
```

## Arguments

- distmat:

  n x n distance matrix

- c:

  Decay parameter

- M:

  n x n projection matrix (I - X(X'X)^-1X')

## Value

Residual-projected covariance matrix
