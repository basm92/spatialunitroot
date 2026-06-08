# Get cbar: find c such that mean(exp(-c \* lower_tri(dist))) = rhobar

Calls the C++ implementation (which works fine).

## Usage

``` r
getcbar_r(rhobar, distmat)
```

## Arguments

- rhobar:

  Target correlation

- distmat:

  n x n distance matrix

## Value

cbar
