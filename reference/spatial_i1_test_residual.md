# Spatial I(1) residual unit root test

Tests residuals from regression Y on X for spatial unit root.

## Usage

``` r
spatial_i1_test_residual(Y, X_in, distmat, emat, q)
```

## Arguments

- Y:

  n x 1 response vector

- X_in:

  n x k design matrix (including intercept)

- distmat:

  n x n normalized distance matrix

- emat:

  q x nrep matrix of standard normal draws

- q:

  Number of low-frequency components

## Value

List with LR, pvalue, ha_parm, cv_vec
