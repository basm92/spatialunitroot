# Package index

## Spatial Unit Root Tests

Diagnostic tests for spatial unit roots

- [`spurtest()`](https://basm92.github.io/spatialunitroot/reference/spurtest.md)
  [`spurtest_i1()`](https://basm92.github.io/spatialunitroot/reference/spurtest.md)
  [`spurtest_i0()`](https://basm92.github.io/spatialunitroot/reference/spurtest.md)
  [`spurtest_i1resid()`](https://basm92.github.io/spatialunitroot/reference/spurtest.md)
  [`spurtest_i0resid()`](https://basm92.github.io/spatialunitroot/reference/spurtest.md)
  : Spatial Unit Root Diagnostic Tests

## Spatial Transformations

Spatial differencing to remove unit roots

- [`spurtransform()`](https://basm92.github.io/spatialunitroot/reference/spurtransform.md)
  : Spatial Differencing Transformations

## Spatial Half-Life

Confidence intervals for spatial persistence

- [`spurhalflife()`](https://basm92.github.io/spatialunitroot/reference/spurhalflife.md)
  : Confidence Intervals for Spatial Half-Life

## Datasets

Example datasets

- [`chetty`](https://basm92.github.io/spatialunitroot/reference/chetty.md)
  : Chetty et al. (2014) Commuting Zone Data

## Package

Package information

- [`spatialunitroot`](https://basm92.github.io/spatialunitroot/reference/spatialunitroot-package.md)
  [`spatialunitroot-package`](https://basm92.github.io/spatialunitroot/reference/spatialunitroot-package.md)
  : spatialunitroot: Spatial Unit Root Diagnostic Tests and
  Transformations

## Internal

Internal functions (C++ backend and computation)

- [`apply_transform()`](https://basm92.github.io/spatialunitroot/reference/apply_transform.md)
  : Apply transformation matrix to a variable
- [`check_missing()`](https://basm92.github.io/spatialunitroot/reference/check_missing.md)
  : Check for missing values in a matrix and report
- [`cluster_matrix()`](https://basm92.github.io/spatialunitroot/reference/cluster_matrix.md)
  : Compute cluster transformation matrix
- [`compute_distances()`](https://basm92.github.io/spatialunitroot/reference/compute_distances.md)
  : Compute distance matrix from coordinates
- [`double_center()`](https://basm92.github.io/spatialunitroot/reference/double_center.md)
  : Double-center a matrix (subtract row means, then column means)
- [`euclidean_distances()`](https://basm92.github.io/spatialunitroot/reference/euclidean_distances.md)
  : Compute Euclidean distance matrix
- [`find_ha_parm_i0()`](https://basm92.github.io/spatialunitroot/reference/find_ha_parm_i0.md)
  : Find ha_parm g for I(0) test
- [`find_ha_parm_i1()`](https://basm92.github.io/spatialunitroot/reference/find_ha_parm_i1.md)
  : Find ha_parm for I(1) test (R implementation)
- [`get_R()`](https://basm92.github.io/spatialunitroot/reference/get_R.md)
  : Get top q eigenvectors of a symmetric matrix
- [`getcbar()`](https://basm92.github.io/spatialunitroot/reference/getcbar.md)
  : Find c such that mean(exp(-c \* lower_tri(distmat))) equals rhobar
- [`getcbar_r()`](https://basm92.github.io/spatialunitroot/reference/getcbar_r.md)
  : Get cbar: find c such that mean(exp(-c \* lower_tri(dist))) = rhobar
- [`getpow_qf()`](https://basm92.github.io/spatialunitroot/reference/getpow_qf.md)
  : Compute power of quadratic form test
- [`getpow_qf_r()`](https://basm92.github.io/spatialunitroot/reference/getpow_qf_r.md)
  : Compute power of quadratic form test (R implementation)
- [`haversine_distances()`](https://basm92.github.io/spatialunitroot/reference/haversine_distances.md)
  : Compute Haversine great-circle distance matrix
- [`iso_matrix()`](https://basm92.github.io/spatialunitroot/reference/iso_matrix.md)
  : Compute isotropic transformation matrix
- [`lbm_gls_matrix()`](https://basm92.github.io/spatialunitroot/reference/lbm_gls_matrix.md)
  : Compute LBM-GLS transformation matrix
- [`lvech()`](https://basm92.github.io/spatialunitroot/reference/lvech.md)
  : Half-vectorization: extract lower triangular elements (excluding
  diagonal)
- [`nn_matrix()`](https://basm92.github.io/spatialunitroot/reference/nn_matrix.md)
  : Compute nearest-neighbor transformation matrix
- [`normalize_distances()`](https://basm92.github.io/spatialunitroot/reference/normalize_distances.md)
  : Normalize distance matrix so maximum distance is 1
- [`resolve_coords()`](https://basm92.github.io/spatialunitroot/reference/resolve_coords.md)
  : Resolve spatial coordinates from various input formats
- [`sigma_dm()`](https://basm92.github.io/spatialunitroot/reference/sigma_dm.md)
  : Compute demeaned spatial covariance with exponential decay
- [`sigma_lbm()`](https://basm92.github.io/spatialunitroot/reference/sigma_lbm.md)
  : Compute LBM covariance matrix from distance matrix
- [`sigma_lbm_dm()`](https://basm92.github.io/spatialunitroot/reference/sigma_lbm_dm.md)
  : Compute demeaned LBM covariance matrix
- [`sigma_residual()`](https://basm92.github.io/spatialunitroot/reference/sigma_residual.md)
  : Compute residual-projected spatial covariance
- [`spatial_i0_core()`](https://basm92.github.io/spatialunitroot/reference/spatial_i0_core.md)
  : Spatial I(0) test (R implementation)
- [`spatial_i0_resid_core()`](https://basm92.github.io/spatialunitroot/reference/spatial_i0_resid_core.md)
  : Spatial I(0) residual test (R implementation)
- [`spatial_i0_test()`](https://basm92.github.io/spatialunitroot/reference/spatial_i0_test.md)
  : Spatial I(0) unit root test
- [`spatial_i0_test_residual()`](https://basm92.github.io/spatialunitroot/reference/spatial_i0_test_residual.md)
  : Spatial I(0) residual unit root test
- [`spatial_i1_core()`](https://basm92.github.io/spatialunitroot/reference/spatial_i1_core.md)
  : Spatial I(1) test (R implementation)
- [`spatial_i1_resid_core()`](https://basm92.github.io/spatialunitroot/reference/spatial_i1_resid_core.md)
  : Spatial I(1) residual test (R implementation)
- [`spatial_i1_test()`](https://basm92.github.io/spatialunitroot/reference/spatial_i1_test.md)
  : Spatial I(1) unit root test
- [`spatial_i1_test_residual()`](https://basm92.github.io/spatialunitroot/reference/spatial_i1_test_residual.md)
  : Spatial I(1) residual unit root test
- [`spatial_persistence()`](https://basm92.github.io/spatialunitroot/reference/spatial_persistence.md)
  : Compute spatial half-life confidence interval
- [`spatial_persistence_core()`](https://basm92.github.io/spatialunitroot/reference/spatial_persistence_core.md)
  : Spatial half-life confidence set (R implementation)
