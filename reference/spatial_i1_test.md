# Spatial I(1) unit root test

Tests H0: spatial unit root (I(1)) vs H1: not I(1).

## Usage

``` r
spatial_i1_test(Y, distmat, emat, q)
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

List with LR, pvalue, ha_parm, cv_vec
