# Compute power of quadratic form test

Compares quadratic forms under null (om0) vs alternative (om1)
covariance. Matches the Stata getpow_qf() function.

## Usage

``` r
getpow_qf(om0_in, om1_in, e)
```

## Arguments

- e:

  q x nrep matrix of standard normal draws

- om0:

  Null hypothesis covariance matrix (q x q)

- om1:

  Alternative hypothesis covariance matrix (q x q)

## Value

Power (scalar between 0 and 1)
