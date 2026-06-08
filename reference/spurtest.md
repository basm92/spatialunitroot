# Spatial Unit Root Diagnostic Tests

## Usage

``` r
spurtest(
  formula,
  data,
  coords = NULL,
  type = c("i1", "i0", "i1resid", "i0resid"),
  latlong = FALSE,
  q = 15,
  nrep = 1e+05,
  seed = NULL,
  ...
)

spurtest_i1(
  formula,
  data,
  coords = NULL,
  latlong = FALSE,
  q = 15,
  nrep = 1e+05,
  seed = NULL,
  ...
)

spurtest_i0(
  formula,
  data,
  coords = NULL,
  latlong = FALSE,
  q = 15,
  nrep = 1e+05,
  seed = NULL,
  ...
)

spurtest_i1resid(
  formula,
  data,
  coords = NULL,
  latlong = FALSE,
  q = 15,
  nrep = 1e+05,
  seed = NULL,
  ...
)

spurtest_i0resid(
  formula,
  data,
  coords = NULL,
  latlong = FALSE,
  q = 15,
  nrep = 1e+05,
  seed = NULL,
  ...
)
```

## Arguments

- formula:

  A formula specifying the variable(s) to test:

  - For simple tests: `y ~ 1` or just the variable name.

  - For residual tests: `y ~ x1 + x2`.

  Alternatively, a
  [`feols`](https://lrberge.github.io/fixest/reference/feols.html) model
  object.

- data:

  A data frame containing the variables.

- coords:

  Coordinate specification: a formula (e.g., `~ s_1 + s_2`), a matrix,
  or `NULL` to auto-detect `s_*` variables.

- type:

  Type of test: `"i1"`, `"i0"`, `"i1resid"`, or `"i0resid"`.

- latlong:

  Logical; are coordinates latitude (s_1) and longitude (s_2) in decimal
  degrees?

- q:

  Number of low-frequency weighted averages. Default: 15.

- nrep:

  Number of Monte Carlo replications. Default: 100000.

- seed:

  Optional seed for reproducibility.

- ...:

  Additional arguments passed to methods.

## Value

An object of class `"spur_test"` (and subclass depending on `type`) with
components:

- statistic:

  LR test statistic.

- p_value:

  Monte Carlo p-value.

- ha_parm:

  Parameter for the alternative hypothesis.

- critical_values:

  Named vector of critical values (1\\ typeType of test performed.
  qNumber of low-frequency components used. nrepNumber of Monte Carlo
  replications. callThe matched call.

Implements four tests from Mueller and Watson (2024) for diagnosing
spatial unit roots: I(1) test, I(0) test, and their residual variants.
Mueller, U. K. and Watson, M. W. (2024). "Spatial Unit Roots and
Spurious Regression." *Econometrica*, 92, 1661–1695.
