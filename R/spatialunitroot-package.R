#' spatialunitroot: Spatial Unit Root Diagnostic Tests and Transformations
#'
#' @useDynLib spatialunitroot, .registration = TRUE
#' @importFrom Rcpp evalCpp
#'
#' R implementation of the spatial unit root diagnostic tests and spatial
#' differencing transformations from Mueller and Watson (2024, Econometrica).
#' Based on the SPUR Stata package by Becker, Boll, and Voth (2025).
#'
#' Provides functions to test for spatial unit roots (I(1) vs I(0) null
#' hypotheses), compute spatial half-life confidence intervals, and apply
#' spatial transformations (LBM-GLS, nearest-neighbor, isotropic, cluster)
#' to remove spatial unit roots.
#'
#' The package integrates with the \pkg{fixest} package for fixed-effects
#' regression.
#'
#' @section Acknowledgements:
#' This R package is based on the SPUR Stata package by Sascha O. Becker,
#' P. David Boll, and Hans-Joachim Voth, who in turn based their code on
#' Matlab replication files from Mueller and Watson (2024).
#'
#' @section Main functions:
#' \describe{
#'   \item{\code{\link{spurtest}}}{Spatial unit root diagnostic tests}
#'   \item{\code{\link{spurtransform}}}{Spatial differencing transformations}
#'   \item{\code{\link{spurhalflife}}}{Half-life confidence intervals}
#' }
#'
#' @references
#' Becker, S. O., Boll, P. D., and Voth, H.-J. (2025).
#' "Spatial Unit Roots in Regressions: A Practitioner's Guide and a Stata
#' Package." \emph{Stata Journal}, forthcoming.
#'
#' Mueller, U. K. and Watson, M. W. (2022).
#' "Spatial Correlation Robust Inference."
#' \emph{Econometrica}, 90, 2901-2935.
#'
#' Mueller, U. K. and Watson, M. W. (2024).
#' "Spatial Unit Roots and Spurious Regression."
#' \emph{Econometrica}, 92, 1661-1695.
#'
#' @author
#' \strong{R implementation & maintainer}: Bas Machielsen \email{bas@machielsen.org}
#'
#' \strong{Original Stata code}: Sascha O. Becker, P. David Boll, Hans-Joachim Voth
#'
#' @docType package
#' @name spatialunitroot-package
NULL

#' Chetty et al. (2014) Commuting Zone Data
#'
#' Data on intergenerational mobility and its correlates across US
#' commuting zones, used in Chetty, Hendren, Kline, and Saez (2014).
#' Adapted from the SPUR Stata package example data.
#'
#' @format A data frame with 707 observations and 39 variables:
#' \describe{
#'   \item{state}{State abbreviation}
#'   \item{cz}{Commuting zone identifier}
#'   \item{s_1}{Latitude (decimal degrees)}
#'   \item{s_2}{Longitude (decimal degrees)}
#'   \item{am}{Absolute upward mobility}
#'   \item{fracblack}{Fraction black}
#'   \item{racseg}{Racial segregation}
#'   \item{segpov25}{Segregation of poverty}
#'   \item{fraccom15}{Fraction commuting < 15 min}
#'   \item{hipc}{Household income per capita}
#'   \item{gini}{Gini coefficient}
#'   \item{incsh1}{Income share top 1\%}
#'   \item{tsr}{Teacher-student ratio}
#'   \item{tsperc}{Test score percentile}
#'   \item{hsdrop}{High school dropout rate}
#'   \item{scind}{Social capital index}
#'   \item{fracrel}{Fraction religious}
#'   \item{crimer}{Crime rate}
#'   \item{fracsm}{Fraction single mothers}
#'   \item{fracdiv}{Fraction divorced}
#'   \item{fracmar}{Fraction married}
#'   \item{loctr}{Local tax rate}
#'   \item{colpc}{College population share}
#'   \item{coltui}{College tuition}
#'   \item{colgrad}{College graduation rate}
#'   \item{manshare}{Manufacturing share}
#'   \item{chimp}{Change in imports}
#'   \item{tlfpr}{Total labor force participation rate}
#'   \item{migirate}{Migration inflow rate}
#'   \item{migorate}{Migration outflow rate}
#'   \item{fracfor}{Fraction foreign-born}
#' }
#'
#' @source
#' Chetty, R., Hendren, N., Kline, P., and Saez, E. (2014).
#' "Where is the Land of Opportunity? The Geography of Intergenerational
#' Mobility in the United States." \emph{Quarterly Journal of Economics},
#' 129(4), 1553-1623.
#'
#' Data sourced from the SPUR Stata package at
#' \url{https://github.com/pdavidboll/SPUR}.
"chetty"
