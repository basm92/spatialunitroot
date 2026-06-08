# Compute LBM-GLS transformation matrix

From sigma_lbm_dm eigendecomposition, keep eigenvalues \> 1e-10, return
V \* diag(1/sqrt(eval)) \* V'.

## Usage

``` r
lbm_gls_matrix(distmat)
```

## Arguments

- distmat:

  n x n distance matrix

## Value

n x n LBM-GLS transformation matrix
