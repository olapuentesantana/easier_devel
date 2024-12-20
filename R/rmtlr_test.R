#' Regularized Multi-Task Linear Regression (RMTLR)
#' model predictions
#'
#' Computes the predictions as a matrix multiplication using
#' both the features input data and the features estimated
#' weights.
#'
#' @importFrom stats na.omit
#'
#' @param x_test numeric matrix containing features values
#' (rows = samples; columns = features).
#' @param coef_matrix numeric matrix containing the parameters
#' values derived from model training (rows = features;
#' columns = tasks).
#'
#' @return Numeric matrix of predicted values (rows = samples;
#' columns = tasks).
#'
#' @examples
#' \donttest{
#' # using a SummarizedExperiment object
#' library(SummarizedExperiment)
#' # Using example exemplary dataset (Mariathasan et al., Nature, 2018)
#' # from easierData. Original processed data is available from
#' # IMvigor210CoreBiologies package.
#' library("easierData")
#'
#' dataset_mariathasan <- easierData::get_Mariathasan2018_PDL1_treatment()
#' RNA_tpm <- assays(dataset_mariathasan)[["tpm"]]
#'
#' # Select a subset of patients to reduce vignette building time.
#' pat_subset <- c(
#'   "SAM76a431ba6ce1", "SAMd3bd67996035", "SAMd3601288319e",
#'   "SAMba1a34b5a060", "SAM18a4dabbc557"
#' )
#' RNA_tpm <- RNA_tpm[, colnames(RNA_tpm) %in% pat_subset]
#'
#' # Computation of TF activity (Garcia-Alonso et al., Genome Res, 2019)
#' tf_activities <- compute_TF_activity(
#'   RNA_tpm = RNA_tpm
#' )
#'
#' # Parameters values should be defined as a matrix
#' # with features as rows and tasks as columns
#' estimated_parameters <- matrix(rnorm(n = (ncol(tf_activities) + 1) * 10),
#'   nrow = ncol(tf_activities) + 1, ncol = 10
#' )
#' rownames(estimated_parameters) <- c("(Intercept)", colnames(tf_activities))
#' colnames(estimated_parameters) <- c(
#'   "CYT", "Ock_IS", "Roh_IS", "chemokines",
#'   "Davoli_IS", "IFNy", "Ayers_expIS", "Tcell_inflamed", "RIR", "TLS"
#' )
#'
#' # Compute predictions using parameters values
#' pred_test <- easier:::rmtlr_test(
#'   x_test = tf_activities,
#'   coef_matrix = estimated_parameters
#' )
#' }
rmtlr_test <- function(x_test,
                       coef_matrix) {

  # Keep intercept
  intercept <- as.matrix(coef_matrix)[1, ]

  # Remove intercept from coef matrix
  coef <- as.matrix(coef_matrix)[-1, , drop = FALSE]

  # match features properly
  pos <- stats::na.omit(match(colnames(x_test), rownames(coef)))
  coef <- coef[pos, , drop = FALSE]

  if (length(coef) > 1) {
    slope <- coef
    fit_pred <- t(matrix(as.matrix(intercept),
      nrow = ncol(slope),
      ncol = nrow(x_test)
    )
    + t(slope) %*% t(as.matrix(x_test[, rownames(slope)])))
  } else {
    slope <- 0
    fit_pred <- matrix(as.matrix(intercept),
      nrow = ncol(slope),
      ncol = nrow(x_test[[1]])
    )
  }
  return(fit_pred)
}
