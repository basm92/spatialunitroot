# spatialunitroot

Spatial Unit Root Diagnostic Tests and Transformations for R

An R implementation of the spatial unit root methods of Müller & Watson
(2024, *Econometrica*), based on the [SPUR Stata
package](https://github.com/pdavidboll/SPUR) by Becker, Boll & Voth
(2025).

**R developer & maintainer**: Bas Machielsen  
**Original Stata code** due to Sascha O. Becker, P. David Boll,
Hans-Joachim Voth

## Installation

``` r
# install.packages("remotes")
remotes::install_github("basm92/spatialunitroot") 
```

## Quick Start

``` r
library(spatialunitroot)
library(fixest)

data(chetty)

# Test for a spatial unit root (I(1) null)
spurtest(am ~ 1, data = chetty, coords = ~ s_1 + s_2,
         type = "i1", latlong = TRUE, seed = 42)
#> Spatial I(1) Test Results
#> ---------------------------------------
#> Test Statistic :  5.759
#> P-value        :  0.388

# Test for NO spatial unit root (I(0) null)
spurtest(am ~ 1, data = chetty, coords = ~ s_1 + s_2,
         type = "i0", latlong = TRUE, seed = 42)
#> Spatial I(0) Test Results
#> ---------------------------------------
#> Test Statistic :  3.044
#> P-value        :  0.001

# Half-life confidence interval
spurhalflife(am ~ 1, data = chetty, coords = ~ s_1 + s_2, latlong = TRUE)
#> Lower bound: 463653 m
#> Upper bound: inf

# Remove spatial unit roots via LBM-GLS transformation
transformed <- spurtransform(~ am + fracblack, data = chetty,
                             coords = ~ s_1 + s_2, prefix = "h_",
                             latlong = TRUE)

# Use with fixest
feols(h_am ~ h_fracblack, data = transformed)
```

## The Three Functions

**[`spurtest()`](https://basm92.github.io/spatialunitroot/reference/spurtest.md)**
— Four spatial unit root tests from Müller & Watson (2024):

| `type` | Null hypothesis | Description |
|----|----|----|
| `"i1"` | Spatial unit root (I(1)) | Tests whether a variable has a spatial unit root |
| `"i0"` | No spatial unit root (I(0)) | Tests whether a variable is spatially stationary |
| `"i1resid"` | Spatial unit root in residuals | Tests regression residuals for a spatial unit root |
| `"i0resid"` | No spatial unit root in residuals | Tests regression residuals for spatial stationarity |

Also accepts `fixest` model objects directly:
`spurtest(feols_model, coords = ~ s_1 + s_2, type = "i1resid")`

**[`spurtransform()`](https://basm92.github.io/spatialunitroot/reference/spurtransform.md)**
— Four spatial differencing transformations: - `"lbmgls"` (default) —
LBM-GLS transformation - `"nn"` — Nearest-neighbor differencing -
`"iso"` — Isotropic differencing (within radius) - `"cluster"` —
Cluster-level demeaning

**[`spurhalflife()`](https://basm92.github.io/spatialunitroot/reference/spurhalflife.md)**
— Confidence sets for the spatial half-life (the distance at which
spatial correlation drops to 1/2).

## Coordinates

Spatial coordinates must be provided as `s_*` variables (for Stata
compatibility) or via a formula:

``` r
coords = ~ s_1 + s_2                 # latitude, longitude (with latlong = TRUE)
coords = ~ x + y                     # Euclidean coordinates
```

## References

- Becker, S. O., Boll, P. D., & Voth, H. J. (2025). Spatial Unit Roots
  in Regressions: A Practitioner’s Guide and a Stata Package. *Stata
  Journal*.
- Müller, U. K. & Watson, M. W. (2024). “Spatial Unit Roots and Spurious
  Regression.” *Econometrica*, 92, 1661–1695.
- Müller, U. K. & Watson, M. W. (2022). “Spatial Correlation Robust
  Inference.” *Econometrica*, 90, 2901–2935.

## License

MIT
