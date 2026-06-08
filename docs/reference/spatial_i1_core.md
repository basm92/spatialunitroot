# Spatial I(1) test (R implementation)

Spatial I(1) test (R implementation)

## Usage

``` r
spatial_i1_core(Y, distmat, emat, q)
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

List with LR, pvalue, ha_parm, cv_vec
