# Resolve spatial coordinates from various input formats

Handles formula, matrix, and auto-detection of s\_\* variables.

## Usage

``` r
resolve_coords(coords, data, latlong = FALSE)
```

## Arguments

- coords:

  A formula (e.g., `~ s_1 + s_2`), a matrix, or NULL for auto-detection
  of s\_\* variables.

- data:

  A data frame containing the coordinate variables.

- latlong:

  Logical; whether coordinates are latitude/longitude.

## Value

A list with components: `coords` (n x d matrix), `latlong` (logical),
`n` (number of observations).
