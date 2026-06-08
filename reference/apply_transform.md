# Apply transformation matrix to a variable

hy = H \* y. For lbmgls, additionally mean-center the result.

## Usage

``` r
apply_transform(y, H, is_lbmgls)
```

## Arguments

- y:

  n x 1 vector

- H:

  n x n transformation matrix

- is_lbmgls:

  Whether to mean-center the result

## Value

n x 1 transformed vector
