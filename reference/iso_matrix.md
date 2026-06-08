# Compute isotropic transformation matrix

iso_mat = I - dist_below_b / rowsum(dist_below_b) where dist_below_b
indicates distances \<= radius

## Usage

``` r
iso_matrix(s, radius, latlong)
```

## Arguments

- s:

  n x d matrix of locations

- radius:

  Radius for differencing

- latlong:

  Whether coordinates are latitude/longitude

## Value

n x n isotropic transformation matrix
