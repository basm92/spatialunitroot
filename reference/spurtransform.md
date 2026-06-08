# Spatial Differencing Transformations

Applies spatial differencing transformations from Mueller and Watson
(2024) to remove spatial unit roots. Supports four transformation types:
LBM-GLS, nearest-neighbor, isotropic, and cluster.

## Usage

``` r
spurtransform(
  formula,
  data,
  coords = NULL,
  prefix = "h_",
  transformation = c("lbmgls", "nn", "iso", "cluster"),
  radius = NULL,
  cluster = NULL,
  latlong = FALSE,
  replace = FALSE,
  separately = FALSE
)
```

## Arguments

- formula:

  A one-sided formula specifying variables to transform, e.g.,
  `~ y + x1 + x2`.

- data:

  A data frame containing the variables.

- coords:

  Coordinate specification: a formula (e.g., `~ s_1 + s_2`), a matrix,
  or `NULL` to auto-detect `s_*` variables.

- prefix:

  Character prefix for transformed variable names. Default: `"h_"`.

- transformation:

  Type of transformation: `"lbmgls"` (default), `"nn"`
  (nearest-neighbor), `"iso"` (isotropic), or `"cluster"`.

- radius:

  Radius for isotropic transformation. Required when
  `transformation = "iso"`.

- cluster:

  Variable name (in `data`) for cluster transformation. Required when
  `transformation = "cluster"`.

- latlong:

  Logical; are coordinates latitude/longitude?

- replace:

  Logical; overwrite existing variables with the prefix? Default:
  `FALSE`.

- separately:

  Logical; apply transformation separately for each variable (handling
  missing values per-variable)?

## Value

The original data frame with added transformed variables (invisibly).

## References

Mueller, U. K. and Watson, M. W. (2024). "Spatial Unit Roots and
Spurious Regression." *Econometrica*, 92, 1661–1695.

## Examples

``` r
if (FALSE) { # \dontrun{
data(chetty)
transformed <- spurtransform(~ am + fracblack, data = chetty,
                              coords = ~ s_1 + s_2, latlong = TRUE,
                              prefix = "h_")
feols(h_am ~ h_fracblack, data = transformed)
} # }
```
