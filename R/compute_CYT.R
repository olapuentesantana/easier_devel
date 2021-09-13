#' Compute cytolytic activity (CYT) score
#'
#' This function calculates the CYT score as the geometric mean of its signature genes.
#'
#' @references Rooney, M.S., Shukla, S.A., Wu, C.J., Getz, G., and Hacohen, N. (2015). Molecular and genetic properties
#' of tumors associated with local immune cytolytic activity. Cell 160, 48–61. https://doi.org/10.1016/j.cell.2014.12.033.
#'
#' @importFrom stats na.omit
#'
#' @param matches numeric vector indicating the index of signature genes in `RNA_tpm`.
#' @param RNA_tpm data.frame containing TPM values with HGNC symbols in rows and samples in columns.
#'
#' @return A numeric matrix with samples in rows and CTY score in a column.
#'
#' @examples
#' # Example does not matter as function will no be exported
compute_CYT <- function(matches, RNA_tpm) {
  # Subset RNA_tpm
  subset_RNA_tpm <- RNA_tpm[matches, ]

  # Calculation: geometric mean (so-called log-average) [TPM, 0.01 offset]
  score <- as.matrix(apply(subset_RNA_tpm + 0.01, 2, function(X) exp(mean(log(X)))))

  return(data.frame(CYT = score, check.names = FALSE))
}
