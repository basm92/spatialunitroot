# Find c such that mean(exp(-c \* lower_tri(distmat))) equals rhobar

Uses exponential bracketing followed by geometric bisection. Matches the
Stata getcbar() function exactly.

## Usage

``` r
getcbar(rhobar, distmat)
```

## Arguments

- rhobar:

  Target correlation

- distmat:

  n x n distance matrix

## Value

The constant cbar
