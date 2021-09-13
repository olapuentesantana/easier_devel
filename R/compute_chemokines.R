#' Compute chemokine signature (chemokines) score
#'
#' Computes chemokines score as the PC1 score that results from applying PCA to z-score expression of its signature genes.
#'
#' @references Messina, J.L., Fenstermacher, D.A., Eschrich, S., Qu, X., Berglund, A.E., Lloyd, M.C., Schell, M.J.,
#' Sondak, V.K., Weber, J.S., and Mulé, J.J. (2012). 12-Chemokine gene signature identifies lymph node-like structures
#' in melanoma: potential for patient selection for immunotherapy? Sci. Rep. 2, 765. https://doi.org/10.1038/srep00765.
#'
#' @importFrom stats na.omit prcomp
#'
#' @param matches numeric vector indicating the index of signature genes in `RNA_tpm`.
#' @param RNA_tpm data.frame containing TPM values with HGNC symbols in rows and samples in columns.
#'
#' @return A numeric matrix with samples in rows and chemokines score in a column.
#'
#' @examples
#'
#' # Example does not matter as function will no be exported
compute_chemokines <- function(matches, RNA_tpm) {
  # Log2 transformation:
  log2_RNA_tpm <- log2(RNA_tpm + 1)

  # Subset gene_expr
  sub_log2_RNA_tpm<- log2_RNA_tpm[matches, ]

  # calculation: using PCA (Z-score calculated within prcomp)
  chemokine_pca <- stats::prcomp(t(sub_log2_RNA_tpm), center = TRUE, scale = TRUE)
  score <- chemokine_pca$x[, 1]

  return(data.frame(chemokines = score, check.names = FALSE))
}
