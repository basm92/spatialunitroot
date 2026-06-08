# Find ha_parm g for I(0) test

Find ha_parm g for I(0) test

## Usage

``` r
find_ha_parm_i0(om_ho, om_i0, om_bm, emat, maxiter = 20)
```

## Arguments

- om_ho:

  q x q null covariance

- om_i0:

  q x q white noise covariance

- om_bm:

  q x q BM covariance

- emat:

  q x nrep standard normal draws

## Value

ha_parm scalar
