# Get top q eigenvectors of a symmetric matrix

Computes eigendecomposition, sorts eigenvalues descending, returns first
q eigenvectors. Matches the Stata get_R() function.

## Usage

``` r
get_R(sigma, q)
```

## Arguments

- sigma:

  n x n symmetric matrix

- q:

  Number of eigenvectors to extract

## Value

n x q matrix of top eigenvectors
