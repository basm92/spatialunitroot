#' @export
print.spur_test <- function(x, ...) {
  cat("\nSpatial ", format_test_type(x$type), " Test Results\n",
      "---------------------------------------\n", sep = "")
  cat("Test Statistic : ", format(x$statistic, digits = 4), "\n")
  cat("P-value        : ", format(x$p_value, digits = 4), "\n")
  cat("---------------------------------------\n\n")

  invisible(x)
}

#' @export
summary.spur_test <- function(object, ...) {
  cat("\nSpatial ", format_test_type(object$type), " Test Summary\n",
      "=======================================\n\n", sep = "")

  cat("Test Statistic : ", format(object$statistic, digits = 6), "\n")
  cat("P-value        : ", format(object$p_value, digits = 6), "\n")
  cat("HA Parameter   : ", format(object$ha_parm, digits = 4), "\n\n")

  cat("Critical Values:\n")
  cv <- object$critical_values
  cat(sprintf("  1%%  : %8.4f\n", cv[1]))
  cat(sprintf("  5%%  : %8.4f\n", cv[2]))
  cat(sprintf("  10%% : %8.4f\n", cv[3]))

  cat("\nTest Configuration:\n")
  cat("  Type        : ", object$type, "\n")
  cat("  q (components) : ", object$q, "\n")
  cat("  MC replications: ", object$nrep, "\n")

  invisible(object)
}

#' @export
plot.spur_test <- function(x, ...) {
  # For I(0) tests, could plot p-values over rho grid
  # For now, provide a simple diagnostic
  cat("Plot method for spur_test is not yet implemented.\n")
  invisible(x)
}

#' @export
print.spur_halflife <- function(x, ...) {
  cat("\nSpatial Half-Life ", x$level * 100, "% Confidence Interval\n",
      "---------------------------------------\n", sep = "")

  if (x$normdist) {
    units <- paste0("fractions of maximum distance ",
                    format(x$max_dist, digits = 4))
  } else {
    units <- "original coordinate units"
  }

  cat("Lower bound: ", format(x$ci_lower, digits = 4), "\n")
  if (is.infinite(x$ci_upper)) {
    cat("Upper bound: inf\n")
  } else {
    cat("Upper bound: ", format(x$ci_upper, digits = 4), "\n")
  }
  cat("---------------------------------------\n")
  cat("(in ", units, ")\n\n", sep = "")

  invisible(x)
}

#' @export
summary.spur_halflife <- function(object, ...) {
  cat("\nSpatial Half-Life Confidence Interval Summary\n",
      "=========================================\n\n", sep = "")

  cat("Confidence level : ", object$level * 100, "%\n")
  cat("CI Lower bound   : ", format(object$ci_lower, digits = 6), "\n")
  if (is.infinite(object$ci_upper)) {
    cat("CI Upper bound   : inf\n")
  } else {
    cat("CI Upper bound   : ", format(object$ci_upper, digits = 6), "\n")
  }
  cat("Max distance     : ", format(object$max_dist, digits = 4), "\n")
  cat("Normalized units : ", object$normdist, "\n")

  cat("\nConfiguration:\n")
  cat("  q (components)  : ", object$q, "\n")
  cat("  MC replications : ", object$nrep, "\n")

  invisible(object)
}

#' @export
plot.spur_halflife <- function(x, ...) {
  cat("Plot method for spur_halflife is not yet implemented.\n")
  invisible(x)
}

#' @export
print.spur_transform <- function(x, ...) {
  cat("\nSpatial Transformation Summary\n",
      "------------------------------\n", sep = "")
  cat("Transformation type:", x$transformation, "\n")
  cat("Variables transformed:", paste(x$variables, collapse = ", "), "\n")
  cat("Prefix:", x$prefix, "\n")
  invisible(x)
}

# Helper: format test type for display
format_test_type <- function(type) {
  switch(type,
    i1       = "I(1)",
    i0       = "I(0)",
    i1resid  = "I(1) Residual",
    i0resid  = "I(0) Residual",
    type
  )
}
