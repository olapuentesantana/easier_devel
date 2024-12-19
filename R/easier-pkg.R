#' easier: predicting immune response
#' by using quantitative descriptors of the
#' tumor microenvironment extracted from RNA-seq
#' data.
#'
#' This package streamlines the assessment of patients'
#' likelihood of immune response using EaSIeR approach.
#'
#' @importFrom stats cor
#' @importFrom rlang .data
#' @importFrom graphics abline lines
#' @importFrom utils View data
#'
#' @references Lapuente-Santana, Oscar, Maisa van Genderen,
#' Peter A. J. Hilbers, Francesca Finotello, and Federica Eduati.
#' 2021. “Interpretable Systems Biomarkers Predict Response to
#' Immune-Checkpoint Inhibitors.” Patterns, 100293.
#' https://doi.org/10.1016/j.patter.2021.100293.
#' 
#' @name easier-pkg
#' @keywords internal
"_PACKAGE"

globalVariables(c("feature", "istop", "threshold", "weight"))
