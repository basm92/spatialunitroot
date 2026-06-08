# Confidence Intervals for Spatial Half-Life

Constructs confidence sets for the half-life of spatial processes based
on Mueller and Watson (2024). The half-life is the distance at which
spatial correlation drops to 1/2.

## Usage

``` r
spurhalflife(
  formula,
  data,
  coords = NULL,
  level = 0.95,
  q = 15,
  nrep = 1e+05,
  latlong = FALSE,
  normdist = FALSE,
  seed = NULL
)
```

## Arguments

- formula:

  A formula specifying the variable, e.g., `y ~ 1`.

- data:

  A data frame containing the variables.

- coords:

  Coordinate specification: a formula (e.g., `~ s_1 + s_2`), a matrix,
  or `NULL` to auto-detect `s_*` variables.

- level:

  Confidence level (between 0 and 1). Default: 0.95.

- q:

  Number of low-frequency weighted averages. Default: 15.

- nrep:

  Number of Monte Carlo replications. Default: 100000.

- latlong:

  Logical; are coordinates latitude/longitude?

- normdist:

  Logical; if `TRUE`, return CI in normalized distance units (fractions
  of maximum distance). If `FALSE` (default), return CI in metres (if
  latlong) or original coordinate units.

- seed:

  Optional seed for reproducibility.

## Value

An object of class `"spur_halflife"` with components:

- ci_lower:

  Lower bound of the confidence interval.

- ci_upper:

  Upper bound of the confidence interval (`Inf` if unbounded).

- max_dist:

  Maximum pairwise distance in original units.

- level:

  Confidence level used.

- normdist:

  Whether CI is in normalized units.

- q:

  Number of low-frequency components.

- nrep:

  Number of Monte Carlo replications.

- call:

  The matched call.

## References

Mueller, U. K. and Watson, M. W. (2024). "Spatial Unit Roots and
Spurious Regression." *Econometrica*, 92, 1661–1695.

## Examples

``` r
if (FALSE) { # \dontrun{
data(chetty)
spurhalflife(am ~ 1, data = chetty, coords = ~ s_1 + s_2, latlong = TRUE)
} # }
```
