#' Resolve spatial coordinates from various input formats
#'
#' Handles formula, matrix, and auto-detection of s_* variables.
#'
#' @param coords A formula (e.g., \code{~ s_1 + s_2}), a matrix, or NULL
#'   for auto-detection of s_* variables.
#' @param data A data frame containing the coordinate variables.
#' @param latlong Logical; whether coordinates are latitude/longitude.
#' @return A list with components: \code{coords} (n x d matrix),
#'   \code{latlong} (logical), \code{n} (number of observations).
#' @keywords internal
resolve_coords <- function(coords, data, latlong = FALSE) {
  if (is.null(coords)) {
    # Auto-detect s_* variables
    s_vars <- grep("^s_[0-9]+$", names(data), value = TRUE)
    if (length(s_vars) == 0) {
      stop("No s_* coordinate variables found. Please specify coords.",
           call. = FALSE)
    }
    s_nums <- as.integer(sub("^s_", "", s_vars))
    s_vars <- s_vars[order(s_nums)]

    # Validate sequential numbering
    for (i in seq_along(s_vars)) {
      expected <- paste0("s_", i)
      if (!expected %in% s_vars) {
        stop("s_* variables not continuously numbered starting from 1",
             call. = FALSE)
      }
    }

    coord_mat <- as.matrix(data[, s_vars, drop = FALSE])
  } else if (inherits(coords, "formula")) {
    # Formula interface: ~ s_1 + s_2
    coord_vars <- all.vars(coords)
    if (length(coord_vars) == 0) {
      stop("coords formula must specify variables, e.g., ~ s_1 + s_2",
           call. = FALSE)
    }
    missing_vars <- setdiff(coord_vars, names(data))
    if (length(missing_vars) > 0) {
      stop("Coordinate variables not found: ",
           paste(missing_vars, collapse = ", "), call. = FALSE)
    }
    coord_mat <- as.matrix(data[, coord_vars, drop = FALSE])
  } else if (is.matrix(coords)) {
    coord_mat <- coords
  } else if (is.data.frame(coords)) {
    coord_mat <- as.matrix(coords)
  } else {
    stop("coords must be a formula, matrix, data frame, or NULL",
         call. = FALSE)
  }

  n <- nrow(coord_mat)
  if (n < 5) {
    stop("Too few locations: need at least 5 observations, found ", n,
         call. = FALSE)
  }

  if (latlong && ncol(coord_mat) != 2) {
    stop("With latlong=TRUE, exactly 2 coordinate columns are required",
         call. = FALSE)
  }

  list(coords = coord_mat, latlong = latlong, n = n)
}

#' Compute distance matrix from coordinates
#'
#' Dispatches to Euclidean or Haversine based on latlong flag.
#'
#' @param coords n x d matrix of coordinates.
#' @param latlong Logical; if TRUE use Haversine great-circle distances.
#' @param normalize Logical; if TRUE normalize so max distance = 1.
#' @return n x n distance matrix.
#' @export
compute_distances <- function(coords, latlong = FALSE, normalize = TRUE) {
  if (latlong) {
    distmat <- haversine_distances(coords)
  } else {
    distmat <- euclidean_distances(coords)
  }
  if (normalize) {
    distmat <- normalize_distances(distmat)
  }
  distmat
}

#' Check for missing values in a matrix and report
#'
#' @param mat A numeric matrix.
#' @param name Character; name for error messages.
#' @keywords internal
check_missing <- function(mat, name = "data") {
  if (anyNA(mat)) {
    n_miss <- sum(is.na(mat))
    stop("Found ", n_miss, " missing value(s) in ", name, call. = FALSE)
  }
  invisible(TRUE)
}
