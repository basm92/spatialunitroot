# spatialunitroot: Spatial Unit Root Diagnostic Tests and Transformations

Implements spatial unit root diagnostic tests and spatial differencing
transformations from Muller and Watson (2024, Econometrica). Provides
functions to test for spatial unit roots (I(1) vs I(0) null hypotheses),
compute spatial half-life confidence intervals, and apply spatial
transformations (LBM-GLS, nearest-neighbor, isotropic, cluster) to
remove spatial unit roots. Integrates with the 'fixest' package for
fixed-effects regression. Based on the SPUR Stata package by Becker,
Boll, and Voth.

## Acknowledgements

This R package is based on the SPUR Stata package by Sascha O. Becker,
P. David Boll, and Hans-Joachim Voth, who in turn based their code on
Matlab replication files from Mueller and Watson (2024).

## Main functions

- [`spurtest`](https://basm92.github.io/spatialunitroot/reference/spurtest.md):

  Spatial unit root diagnostic tests

- [`spurtransform`](https://basm92.github.io/spatialunitroot/reference/spurtransform.md):

  Spatial differencing transformations

- [`spurhalflife`](https://basm92.github.io/spatialunitroot/reference/spurhalflife.md):

  Half-life confidence intervals

## References

Becker, S. O., Boll, P. D., and Voth, H.-J. (2025). "Spatial Unit Roots
in Regressions: A Practitioner's Guide and a Stata Package." *Stata
Journal*, forthcoming.

Mueller, U. K. and Watson, M. W. (2022). "Spatial Correlation Robust
Inference." *Econometrica*, 90, 2901-2935.

Mueller, U. K. and Watson, M. W. (2024). "Spatial Unit Roots and
Spurious Regression." *Econometrica*, 92, 1661-1695.

## See also

Useful links:

- <https://github.com/basm92/spatialunitroot>

- Report bugs at <https://github.com/basm92/spatialunitroot/issues>

## Author

**R implementation & maintainer**: Bas Machielsen <bas@machielsen.org>

**Original Stata code**: Sascha O. Becker, P. David Boll, Hans-Joachim
Voth
