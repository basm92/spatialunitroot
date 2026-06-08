# Find ha_parm for I(1) test (R implementation)

Find ha_parm for I(1) test (R implementation)

## Usage

``` r
find_ha_parm_i1(om_ho, distmat, Rmat, emat, maxiter = 20)
```

## Arguments

- om_ho:

  q x q null covariance

- distmat:

  n x n distance matrix

- Rmat:

  n x q eigenvector matrix

- emat:

  q x nrep standard normal draws

## Value

ha_parm scalar
