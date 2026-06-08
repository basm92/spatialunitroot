# Spatial I(0) test (R implementation)

Spatial I(0) test (R implementation)

## Usage

``` r
spatial_i0_core(Y, distmat, emat, q)
```

## Arguments

- Y:

  n x 1 vector

- distmat:

  n x n normalized distance matrix

- emat:

  q x nrep standard normal draws

- q:

  Number of low-frequency components

## Value

List with LR, pvalue, cvalue, ha_parm, rho_grid, pvalue_mat, cvalue_mat
