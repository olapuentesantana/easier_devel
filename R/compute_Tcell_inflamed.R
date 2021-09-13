#' Compute T cell-inflamed signature (Tcell_inflamed) score
#'
#' This function calculates Tcell_inflamed score as a weighted sum of housekeeping normalized expression of its signature genes.
#' Weightes were available at Table S2B from Cristescu R, et al. Pan-tumor genomic biomarkers for PD-1 checkpoint
#' blockade-based immunotherapy. Science. (2018) 362:eaar3593. doi: 10.1126/science.aar3593.
#'
#' @references Ayers, M., Lunceford, J., Nebozhyn, M., Murphy, E., Loboda, A., Kaufman, D.R., Albright,
#' A., Cheng, J.D., Kang, S.P., Shankaran, V., et al. (2017). IFN-γ-related mRNA profile predicts clinical
#' response to PD-1 blockade. J. Clin. Invest. 127, 2930–2940. https://doi.org/10.1172/JCI91190.
#'
#' @importFrom stats na.omit
#'
#' @param housekeeping numeric vector indicating the index of houskeeping genes in `RNA_tpm`.
#' @param predictors numeric vector indicating the index of predictor genes in `RNA_tpm`.
#' @param weights numeric vector containing the weights.
#' @param RNA_tpm data.frame containing TPM values with HGNC symbols in rows and samples in columns.
#'
#' @return A numeric matrix with samples in rows and Tcell_inflamed score in a column.
#'
#' @examples
#'
#' # Example does not matter as function will no be exported
compute_Tcell_inflamed <- function(housekeeping, predictors, weights, RNA_tpm) {

  # Log2 transformation:
  log2_RNA_tpm <- log2(RNA_tpm + 1)

  # Subset log2.RNA_tpm
  ## housekeeping
  log2_RNA_tpm_housekeeping <- log2_RNA_tpm[housekeeping, ]
  ## predictors
  log2_RNA_tpm_predictors <- log2_RNA_tpm[predictors, ]
  weights <- weights[rownames(log2_RNA_tpm_predictors)]

  # Housekeeping normalization
  average_log2_RNA_tpm_housekeeping <- apply(log2_RNA_tpm_housekeeping, 2, mean)
  log2_RNA_tpm_predictors_norm <- sweep(log2_RNA_tpm_predictors, 2, average_log2_RNA_tpm_housekeeping, FUN = "-")

  # Calculation: weighted sum of the normalized predictor gene values
  tidy <- match(rownames(log2_RNA_tpm_predictors_norm), names(weights))

  # Transform vector to matrix
  weights <- matrix(weights, ncol = 1, dimnames = list(names(weights)))
  score <- t(log2_RNA_tpm_predictors_norm[tidy, ]) %*% weights

  return(data.frame(Tcell_inflamed = score, check.names = FALSE))
}
