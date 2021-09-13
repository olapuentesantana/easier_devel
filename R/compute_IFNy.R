#' Compute IFNy signature (IFNy) score
#'
#' This function calculates IFNy signature score as the average expression of its signature genes.
#'
#' @references Ayers, M., Lunceford, J., Nebozhyn, M., Murphy, E., Loboda, A., Kaufman, D.R., Albright,
#' A., Cheng, J.D., Kang, S.P., Shankaran, V., et al. (2017). IFN-γ-related mRNA profile predicts clinical
#' response to PD-1 blockade. J. Clin. Invest. 127, 2930–2940. https://doi.org/10.1172/JCI91190.
#'
#' @importFrom stats na.omit
#'
#' @param matches numeric vector indicating the index of signature genes in `RNA_tpm`.
#' @param RNA_tpm data.frame containing TPM values with HGNC symbols in rows and samples in columns.
#'
#' @return A numeric matrix with samples in rows and IFNy score in a column.
#'
#' @examples
#'
#' # Example does not matter as function will no be exported
compute_IFNy <- function(matches, RNA_tpm) {
  # Log2 transformation:
  log2_RNA_tpm <- log2(RNA_tpm + 1)

  # Subset log2.RNA_tpm
  sub_log2_RNA_tpm <- log2_RNA_tpm[matches, ]

  # Calculation: average of the included genes for the IFN-y signature
  score <- apply(sub_log2_RNA_tpm, 2, mean)

  return(data.frame(IFNy = score, check.names = FALSE))
}
