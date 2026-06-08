#' Spatial Differencing Transformations
#'
#' Applies spatial differencing transformations from Mueller and Watson (2024)
#' to remove spatial unit roots. Supports four transformation types:
#' LBM-GLS, nearest-neighbor, isotropic, and cluster.
#'
#' @param formula A one-sided formula specifying variables to transform,
#'   e.g., \code{~ y + x1 + x2}.
#' @param data A data frame containing the variables.
#' @param coords Coordinate specification: a formula (e.g.,
#'   \code{~ s_1 + s_2}), a matrix, or \code{NULL} to auto-detect
#'   \code{s_*} variables.
#' @param prefix Character prefix for transformed variable names.
#'   Default: \code{"h_"}.
#' @param transformation Type of transformation: \code{"lbmgls"} (default),
#'   \code{"nn"} (nearest-neighbor), \code{"iso"} (isotropic), or
#'   \code{"cluster"}.
#' @param radius Radius for isotropic transformation. Required when
#'   \code{transformation = "iso"}.
#' @param cluster Variable name (in \code{data}) for cluster transformation.
#'   Required when \code{transformation = "cluster"}.
#' @param latlong Logical; are coordinates latitude/longitude?
#' @param replace Logical; overwrite existing variables with the prefix?
#'   Default: \code{FALSE}.
#' @param separately Logical; apply transformation separately for each
#'   variable (handling missing values per-variable)?
#'
#' @return The original data frame with added transformed variables
#'   (invisibly).
#'
#' @references
#' Mueller, U. K. and Watson, M. W. (2024).
#' "Spatial Unit Roots and Spurious Regression."
#' \emph{Econometrica}, 92, 1661--1695.
#'
#' @examples
#' \dontrun{
#' data(chetty)
#' transformed <- spurtransform(~ am + fracblack, data = chetty,
#'                               coords = ~ s_1 + s_2, latlong = TRUE,
#'                               prefix = "h_")
#' feols(h_am ~ h_fracblack, data = transformed)
#' }
#'
#' @export
spurtransform <- function(formula, data, coords = NULL,
                          prefix = "h_",
                          transformation = c("lbmgls", "nn", "iso", "cluster"),
                          radius = NULL, cluster = NULL,
                          latlong = FALSE, replace = FALSE,
                          separately = FALSE) {

  transformation <- match.arg(transformation)
  cl <- match.call()

  # Validate transformation options
  if (transformation == "iso" && is.null(radius)) {
    stop("Radius required for isotropic transformation", call. = FALSE)
  }
  if (transformation == "iso" && radius <= 0) {
    stop("Radius must be positive", call. = FALSE)
  }
  if (transformation == "cluster" && is.null(cluster)) {
    stop("Cluster variable required for cluster transformation",
         call. = FALSE)
  }
  if (transformation != "iso" && !is.null(radius)) {
    stop("Radius only allowed with transformation='iso'", call. = FALSE)
  }
  if (transformation != "cluster" && !is.null(cluster)) {
    stop("Cluster only allowed with transformation='cluster'", call. = FALSE)
  }

  # Parse variables from formula
  if (!inherits(formula, "formula")) {
    stop("formula must be a one-sided formula, e.g., ~ y + x1 + x2",
         call. = FALSE)
  }
  vars <- all.vars(formula)
  if (length(vars) == 0) {
    stop("No variables specified in formula", call. = FALSE)
  }

  # Check variables exist
  missing_vars <- setdiff(vars, names(data))
  if (length(missing_vars) > 0) {
    stop("Variables not found: ", paste(missing_vars, collapse = ", "),
         call. = FALSE)
  }

  # Resolve coordinates
  coord_info <- resolve_coords(coords, data, latlong)
  coord_mat <- coord_info$coords

  # Handle cluster variable
  if (transformation == "cluster") {
    cluster_vec <- data[[cluster]]
    if (is.null(cluster_vec)) {
      stop("Cluster variable '", cluster, "' not found", call. = FALSE)
    }
    # Convert to integer cluster IDs
    cluster_vec <- as.integer(as.factor(cluster_vec))
  }

  # Process each variable
  for (v in vars) {
    new_name <- paste0(prefix, v)

    # Check if new variable already exists
    if (new_name %in% names(data) && !replace) {
      stop("Variable '", new_name, "' already exists. Use replace=TRUE ",
           "to overwrite.", call. = FALSE)
    }

    if (separately) {
      # Find non-missing observations for this variable
      ok <- complete.cases(data[[v]], coord_mat)
      if (transformation == "cluster") {
        ok <- ok & !is.na(cluster_vec)
      }
    } else {
      # Use all complete cases across all variables and coordinates
      all_data <- as.matrix(data[, vars, drop = FALSE])
      if (transformation == "cluster") {
        ok <- complete.cases(all_data, coord_mat, cluster_vec)
      } else {
        ok <- complete.cases(all_data, coord_mat)
      }
    }

    if (sum(ok) < 5) {
      warning("Variable '", v, "': too few complete observations, skipping")
      next
    }

    y <- data[[v]][ok]
    cmat <- coord_mat[ok, , drop = FALSE]

    # Build transformation matrix
    if (transformation == "lbmgls") {
      distmat <- compute_distances(cmat, latlong = latlong, normalize = TRUE)
      H <- lbm_gls_matrix(distmat)
      hy <- apply_transform(as.vector(y), H, TRUE)
    } else if (transformation == "nn") {
      H <- nn_matrix(cmat, latlong)
      hy <- apply_transform(as.vector(y), H, FALSE)
    } else if (transformation == "iso") {
      H <- iso_matrix(cmat, radius, latlong)
      hy <- apply_transform(as.vector(y), H, FALSE)
    } else if (transformation == "cluster") {
      cl_vec <- cluster_vec[ok]
      H <- cluster_matrix(as.vector(cl_vec))
      hy <- apply_transform(as.vector(y), H, FALSE)
    }

    # Store result
    data[[new_name]] <- NA_real_
    data[[new_name]][ok] <- as.vector(hy)
  }

  invisible(data)
}
