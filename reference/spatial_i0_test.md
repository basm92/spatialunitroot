# Spatial I(0) unit root test

Tests H0: no spatial unit root (I(0)) vs H1: spatial unit root. Uses a
grid of rho values for the null hypothesis.

## Usage

``` r
spatial_i0_test(Y, distmat, emat, q)
```

## Arguments

- Y:

  n x 1 vector of observations

- distmat:

  n x n normalized distance matrix

- emat:

  q x nrep matrix of standard normal draws

- q:

  Number of low-frequency components

## Value

List with LR, pvalue, cvalue, ha_parm, rho_grid, pvalue_mat, cvalue_mat
