# Compute distance matrix from coordinates

Dispatches to Euclidean or Haversine based on latlong flag.

## Usage

``` r
compute_distances(coords, latlong = FALSE, normalize = TRUE)
```

## Arguments

- coords:

  n x d matrix of coordinates.

- latlong:

  Logical; if TRUE use Haversine great-circle distances.

- normalize:

  Logical; if TRUE normalize so max distance = 1.

## Value

n x n distance matrix.
