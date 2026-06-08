# Spatial I(0) residual test (R implementation)

Spatial I(0) residual test (R implementation)

## Usage

``` r
spatial_i0_resid_core(Y, Xmat, distmat, emat, q)
```

## Arguments

- Y:

  n x 1 response vector

- Xmat:

  n x k design matrix (including intercept)

- distmat:

  n x n normalized distance matrix

- emat:

  q x nrep standard normal draws

- q:

  Number of low-frequency components

## Value

List with LR, pvalue, cvalue, ha_parm, rho_grid, pvalue_mat, cvalue_mat
