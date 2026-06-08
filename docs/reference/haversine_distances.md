# Compute Haversine great-circle distance matrix

Distances are divided by pi to match Stata convention (fraction of
semi-circle, so max distance = 1).

## Usage

``` r
haversine_distances(coords)
```

## Arguments

- coords:

  n x 2 matrix with latitude in column 0, longitude in column 1

## Value

n x n Haversine distance matrix (scaled by 1/pi)
