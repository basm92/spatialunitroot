#' Spatial Unit Root Diagnostic Tests
#'
#' Implements four tests from Mueller and Watson (2024) for diagnosing
#' spatial unit roots: I(1) test, I(0) test, and their residual variants.
#'
#' @param formula A formula specifying the variable(s) to test:
#'   \itemize{
#'     \item For simple tests: \code{y ~ 1} or just the variable name.
#'     \item For residual tests: \code{y ~ x1 + x2}.
#'   }
#'   Alternatively, a \code{\link[fixest]{feols}} model object.
#' @param data A data frame containing the variables.
#' @param coords Coordinate specification: a formula (e.g.,
#'   \code{~ s_1 + s_2}), a matrix, or \code{NULL} to auto-detect
#'   \code{s_*} variables.
#' @param type Type of test: \code{"i1"}, \code{"i0"}, \code{"i1resid"},
#'   or \code{"i0resid"}.
#' @param latlong Logical; are coordinates latitude (s_1) and longitude
#'   (s_2) in decimal degrees?
#' @param q Number of low-frequency weighted averages. Default: 15.
#' @param nrep Number of Monte Carlo replications. Default: 100000.
#' @param seed Optional seed for reproducibility.
#' @param ... Additional arguments passed to methods.
#'
#' @return An object of class \code{"spur_test"} (and subclass depending on
#'   \code{type}) with components:
#'   \item{statistic}{LR test statistic.}
#'   \item{p_value}{Monte Carlo p-value.}
#'   \item{ha_parm}{Parameter for the alternative hypothesis.}
#'   \item{critical_values}{Named vector of critical values (1\%, 5\%, 10\%).}
#'   \item{type}{Type of test performed.}
#'   \item{q}{Number of low-frequency components used.}
#'   \item{nrep}{Number of Monte Carlo replications.}
#'   \item{call}{The matched call.}
#'
#' @references
#' Mueller, U. K. and Watson, M. W. (2024).
#' "Spatial Unit Roots and Spurious Regression."
#' \emph{Econometrica}, 92, 1661--1695.
#'
#' @examples
#' \dontrun{
#' data(chetty)
#' spurtest(am ~ 1, data = chetty, coords = ~ s_1 + s_2,
#'          type = "i1", latlong = TRUE)
#' }
#'
#' @export
spurtest <- function(formula, data, coords = NULL,
                     type = c("i1", "i0", "i1resid", "i0resid"),
                     latlong = FALSE, q = 15, nrep = 100000,
                     seed = NULL, ...) {
  UseMethod("spurtest")
}

#' @export
#' @method spurtest default
spurtest.default <- function(formula, data, coords = NULL,
                              type = c("i1", "i0", "i1resid", "i0resid"),
                              latlong = FALSE, q = 15, nrep = 100000,
                              seed = NULL, ...) {
  type <- match.arg(type)
  cl <- match.call()

  # Parse formula
  if (!inherits(formula, "formula")) {
    formula <- as.formula(paste(deparse(substitute(formula)), "~ 1"))
  }

  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Determine if this is a residual test
  has_covariates <- type %in% c("i1resid", "i0resid")

  # Build model frame
  mf <- model.frame(formula, data = data, na.action = na.pass)
  y <- model.response(mf, "numeric")
  if (is.null(y)) stop("Response variable not found in formula", call. = FALSE)

  # Get covariate matrix for residual tests
  if (has_covariates) {
    mt <- terms(formula, data = data)
    x_vars <- attr(mt, "term.labels")
    if (length(x_vars) == 0) {
      # No covariates beyond intercept -> use simple test
      has_covariates <- FALSE
      if (type == "i1resid") type <- "i1" else type <- "i0"
    } else {
      mf_x <- model.frame(formula, data = data, na.action = na.pass)
      xmat <- model.matrix(formula, data = mf_x)
    }
  }

  # Resolve coordinates
  coord_info <- resolve_coords(coords, data, latlong)
  coord_mat <- coord_info$coords

  # Find complete cases
  if (has_covariates) {
    ok <- complete.cases(y, coord_mat, xmat)
    y <- y[ok]
    coord_mat <- coord_mat[ok, , drop = FALSE]
    xmat <- xmat[ok, , drop = FALSE]
  } else {
    ok <- complete.cases(y, coord_mat)
    y <- y[ok]
    coord_mat <- coord_mat[ok, , drop = FALSE]
  }

  if (length(y) < 5) {
    stop("Too few complete observations: need at least 5, found ",
         length(y), call. = FALSE)
  }

  # Check for missing values in coordinates
  if (anyNA(coord_mat)) {
    stop("Missing values in coordinate variables", call. = FALSE)
  }

  # Compute distance matrix
  distmat <- compute_distances(coord_mat, latlong = latlong, normalize = TRUE)

  # Generate Monte Carlo draws
  emat <- matrix(stats::rnorm(q * nrep), nrow = q, ncol = nrep)

  # Run test via R implementation (more robust than C++ for eigendecomposition)
  if (type == "i1") {
    result <- spatial_i1_core(as.vector(y), distmat, emat, q)
  } else if (type == "i0") {
    result <- spatial_i0_core(as.vector(y), distmat, emat, q)
  } else if (type == "i1resid") {
    result <- spatial_i1_resid_core(as.vector(y), xmat, distmat, emat, q)
  } else if (type == "i0resid") {
    result <- spatial_i0_resid_core(as.vector(y), xmat, distmat, emat, q)
  }

  # Build S3 object
  cv_vals <- if (!is.null(result$cv)) result$cv else result$cvalue
  structure(
    list(
      statistic = result$LR,
      p_value   = result$pvalue,
      ha_parm   = result$ha_parm,
      critical_values = stats::setNames(cv_vals,
                                      c("1%", "5%", "10%")),
      type      = type,
      q         = q,
      nrep      = nrep,
      call      = cl
    ),
    class = c(paste0("spur_test_", type), "spur_test")
  )
}

#' @export
#' @method spurtest fixest
spurtest.fixest <- function(formula, data, coords = NULL,
                             type = c("i1resid", "i0resid"),
                             latlong = FALSE, q = 15, nrep = 100000,
                             seed = NULL, ...) {
  type <- match.arg(type)
  model <- formula
  cl <- match.call()

  # Extract response and model matrix from fixest model
  y <- as.vector(stats::residuals(model))
  xmat <- stats::model.matrix(model, type = "rhs")

  # Get the original data and identify which observations were used
  data_env <- model$call$data
  if (is.name(data_env)) {
    orig_data <- eval(data_env, parent.frame())
  } else if (is.data.frame(data_env)) {
    orig_data <- data_env
  } else {
    stop("Cannot determine the data used in the fixest model.", call. = FALSE)
  }

  # Get observation indices: find which rows are complete for the model variables
  all_vars <- all.vars(model$fml)
  ok_orig <- complete.cases(orig_data[, intersect(all_vars, names(orig_data)), drop = FALSE])
  obs_used <- which(ok_orig)

  # Get coordinates
  if (is.null(coords)) {
    coord_info <- resolve_coords(NULL, orig_data, latlong)
  } else if (inherits(coords, "formula")) {
    coord_info <- resolve_coords(coords, orig_data, latlong)
  } else if (is.matrix(coords) || is.data.frame(coords)) {
    coord_info <- list(coords = as.matrix(coords), latlong = latlong, n = nrow(coords))
  } else {
    coord_info <- resolve_coords(coords, orig_data, latlong)
  }

  # Subset coordinates to match model observations
  coord_mat <- coord_info$coords
  if (length(obs_used) <= nrow(coord_mat) && length(obs_used) == length(y)) {
    coord_mat <- coord_mat[obs_used, , drop = FALSE]
  }

  # Remove any remaining NAs in coordinates
  ok <- complete.cases(coord_mat)
  y <- y[ok]
  coord_mat <- coord_mat[ok, , drop = FALSE]
  xmat <- xmat[ok, , drop = FALSE]

  if (length(y) < 5) {
    stop("Too few complete observations: need at least 5", call. = FALSE)
  }

  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Compute distance matrix
  distmat <- compute_distances(coord_mat, latlong = latlong, normalize = TRUE)

  # Generate Monte Carlo draws
  emat <- matrix(stats::rnorm(q * nrep), nrow = q, ncol = nrep)

  # Run test (R implementation)
  if (type == "i1resid") {
    result <- spatial_i1_resid_core(as.vector(y), xmat, distmat, emat, q)
  } else {
    result <- spatial_i0_resid_core(as.vector(y), xmat, distmat, emat, q)
  }

  # Build S3 object
  cv_vals <- if (!is.null(result$cv)) result$cv else result$cvalue
  structure(
    list(
      statistic = result$LR,
      p_value   = result$pvalue,
      ha_parm   = result$ha_parm,
      critical_values = stats::setNames(cv_vals,
                                      c("1%", "5%", "10%")),
      type      = type,
      q         = q,
      nrep      = nrep,
      call      = cl
    ),
    class = c(paste0("spur_test_", type), "spur_test")
  )
}

#' @rdname spurtest
#' @export
spurtest_i1 <- function(formula, data, coords = NULL, latlong = FALSE,
                        q = 15, nrep = 100000, seed = NULL, ...) {
  spurtest(formula, data = data, coords = coords, type = "i1",
           latlong = latlong, q = q, nrep = nrep, seed = seed, ...)
}

#' @rdname spurtest
#' @export
spurtest_i0 <- function(formula, data, coords = NULL, latlong = FALSE,
                        q = 15, nrep = 100000, seed = NULL, ...) {
  spurtest(formula, data = data, coords = coords, type = "i0",
           latlong = latlong, q = q, nrep = nrep, seed = seed, ...)
}

#' @rdname spurtest
#' @export
spurtest_i1resid <- function(formula, data, coords = NULL, latlong = FALSE,
                             q = 15, nrep = 100000, seed = NULL, ...) {
  spurtest(formula, data = data, coords = coords, type = "i1resid",
           latlong = latlong, q = q, nrep = nrep, seed = seed, ...)
}

#' @rdname spurtest
#' @export
spurtest_i0resid <- function(formula, data, coords = NULL, latlong = FALSE,
                             q = 15, nrep = 100000, seed = NULL, ...) {
  spurtest(formula, data = data, coords = coords, type = "i0resid",
           latlong = latlong, q = q, nrep = nrep, seed = seed, ...)
}
