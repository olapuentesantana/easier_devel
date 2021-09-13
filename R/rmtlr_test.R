#' Regularized Multi-Task Linear Regression (RMTLR) model predictions
#'
#' Computes the predictions from the results of RMTLR models optimization
#'
#' @importFrom stats na.omit
#'
#' @export
#'
#' @param x_test numeric matrix containing features values (rows = samples; columns = features).
#' @param coef_matrix numeric matrix containing the parameters values derived from model training
#' (rows = features; columns = tasks).
#'
#' @return Numeric matrix of predicted values (rows = samples; columns = tasks).
#'
#' @examples
#' # Load exemplary dataset (Mariathasan et al., Nature, 2018) from easierData.
#' # Original processed data is available from IMvigor210CoreBiologies package.
#' library("easierData")
#' dataset_mariathasan <- easierData::get_Mariathasan2018_PDL1_treatment()
#' RNA_tpm <- dataset_mariathasan@assays@data@listData[["tpm"]]
#'
#' # Computation of cell fractions (Finotello et al., Genome Med, 2019)
#' cell_fractions <- compute_cell_fractions(RNA_tpm = RNA_tpm)
#'
#' # Parameters values should be defined with features as rows and tasks as columns
#' estimated_parameters <- matrix(rnorm(n = (ncol(cell_fractions) + 1) * 10),
#'   nrow = ncol(cell_fractions) + 1, ncol = 10
#' )
#' rownames(estimated_parameters) <- c("(Intercept)", colnames(cell_fractions))
#' colnames(estimated_parameters) <- c(
#'   "CYT", "Ock_IS", "Roh_IS", "chemokines",
#'   "Davoli_IS", "IFNy", "Ayers_expIS", "Tcell_inflamed", "RIR", "TLS"
#' )
#'
#' # Compute predictions using parameters values
#' pred_test <- rmtlr_test(
#'   x_test = cell_fractions,
#'   coef_matrix = estimated_parameters
#' )
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
    fit_pred <- t(matrix(as.matrix(intercept), nrow = ncol(slope), ncol = nrow(x_test))
    + t(slope) %*% t(as.matrix(x_test[, rownames(slope)])))
  } else {
    slope <- 0
    fit_pred <- matrix(as.matrix(intercept), nrow = ncol(slope), ncol = nrow(x_test[[1]]))
  }
  return(fit_pred)
}
