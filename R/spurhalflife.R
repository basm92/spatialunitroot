#' Confidence Intervals for Spatial Half-Life
#'
#' Constructs confidence sets for the half-life of spatial processes
#' based on Mueller and Watson (2024). The half-life is the distance at
#' which spatial correlation drops to 1/2.
#'
#' @param formula A formula specifying the variable, e.g., \code{y ~ 1}.
#' @param data A data frame containing the variables.
#' @param coords Coordinate specification: a formula (e.g.,
#'   \code{~ s_1 + s_2}), a matrix, or \code{NULL} to auto-detect
#'   \code{s_*} variables.
#' @param level Confidence level (between 0 and 1). Default: 0.95.
#' @param q Number of low-frequency weighted averages. Default: 15.
#' @param nrep Number of Monte Carlo replications. Default: 100000.
#' @param latlong Logical; are coordinates latitude/longitude?
#' @param normdist Logical; if \code{TRUE}, return CI in normalized
#'   distance units (fractions of maximum distance). If \code{FALSE}
#'   (default), return CI in metres (if latlong) or original coordinate units.
#' @param seed Optional seed for reproducibility.
#'
#' @return An object of class \code{"spur_halflife"} with components:
#'   \item{ci_lower}{Lower bound of the confidence interval.}
#'   \item{ci_upper}{Upper bound of the confidence interval
#'     (\code{Inf} if unbounded).}
#'   \item{max_dist}{Maximum pairwise distance in original units.}
#'   \item{level}{Confidence level used.}
#'   \item{normdist}{Whether CI is in normalized units.}
#'   \item{q}{Number of low-frequency components.}
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
#' spurhalflife(am ~ 1, data = chetty, coords = ~ s_1 + s_2, latlong = TRUE)
#' }
#'
#' @export
spurhalflife <- function(formula, data, coords = NULL,
                         level = 0.95, q = 15, nrep = 100000,
                         latlong = FALSE, normdist = FALSE,
                         seed = NULL) {
  cl <- match.call()

  if (level >= 1 || level <= 0) {
    stop("level must be between 0 and 1", call. = FALSE)
  }

  # Parse formula
  if (!inherits(formula, "formula")) {
    formula <- as.formula(paste(deparse(substitute(formula)), "~ 1"))
  }

  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Get response variable
  mf <- model.frame(formula, data = data, na.action = na.pass)
  y <- model.response(mf, "numeric")
  if (is.null(y)) stop("Response variable not found in formula", call. = FALSE)

  # Resolve coordinates
  coord_info <- resolve_coords(coords, data, latlong)
  coord_mat <- coord_info$coords

  # Complete cases
  ok <- complete.cases(y, coord_mat)
  y <- y[ok]
  coord_mat <- coord_mat[ok, , drop = FALSE]

  if (length(y) < 5) {
    stop("Too few complete observations: need at least 5, found ",
         length(y), call. = FALSE)
  }

  # Compute distance matrix (normalized for internal computation)
  distmat_norm <- compute_distances(coord_mat, latlong = latlong,
                                    normalize = TRUE)

  # Generate Monte Carlo draws
  emat <- matrix(stats::rnorm(q * nrep), nrow = q, ncol = nrep)

  # Run persistence computation (R implementation)
  result <- spatial_persistence_core(as.vector(y), distmat_norm, emat, level)

  ci_lower <- result$ci_lower
  ci_upper <- result$ci_upper

  # Compute max distance in original units
  distmat_raw <- compute_distances(coord_mat, latlong = latlong,
                                    normalize = FALSE)
  if (latlong) {
    # Great circle distance in metres
    max_dist <- max(distmat_raw) * pi * 6371000.009 * 2
  } else {
    max_dist <- max(distmat_raw)
  }

  # Convert CI from normalized to original units if not normdist
  if (!normdist && !is.na(ci_lower) && !is.na(ci_upper)) {
    ci_lower <- ci_lower * max_dist
    ci_upper <- ci_upper * max_dist
  }

  # Handle unbounded upper CI (>= 100 in normalized units)
  if (!is.na(ci_upper) && ci_upper >= 100 * max_dist && !normdist) {
    ci_upper <- Inf
  }
  if (!is.na(ci_upper) && ci_upper >= 100 && normdist) {
    ci_upper <- Inf
  }

  structure(
    list(
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      max_dist = max_dist,
      level    = level,
      normdist = normdist,
      q        = q,
      nrep     = nrep,
      call     = cl
    ),
    class = "spur_halflife"
  )
}
