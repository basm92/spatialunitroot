# Compute cluster transformation matrix

clust_mat = I - (clust_id == clust_id') / rowsum(clust_id == clust_id')

## Usage

``` r
cluster_matrix(cluster)
```

## Arguments

- cluster:

  n x 1 vector of integer cluster IDs

## Value

n x n cluster transformation matrix
