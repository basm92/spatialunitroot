# Compute nearest-neighbor transformation matrix

NN_mat = I - rowmins_mat, where each row has 1/k for the k nearest
neighbors.

## Usage

``` r
nn_matrix(s, latlong)
```

## Arguments

- s:

  n x d matrix of locations

- latlong:

  Whether coordinates are latitude/longitude

## Value

n x n NN transformation matrix
