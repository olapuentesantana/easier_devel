#' Compute Roh immune score (Roh_IS)
#'
#' This function computes Roh_IS score as the geometric-mean of its signature genes.
#'
#' @references Roh, W., Chen, P.-L., Reuben, A., Spencer, C.N., Prieto, P.A., Miller, J.P., Gopalakrishnan, V.,
#' Wang, F., Cooper, Z.A., Reddy, S.M., et al. (2017). Integrated molecular analysis of tumor biopsies on sequential
#' CTLA-4 and PD-1 blockade reveals markers of response and resistance. Sci. Transl. Med. 9.
#' https://doi.org/10.1126/scitranslmed.aah3560.
#'
#' @importFrom stats na.omit
#'
#' @param matches numeric vector indicating the index of signature genes in `RNA_tpm`.
#' @param RNA_tpm data.frame containing TPM values with HGNC symbols in rows and samples in columns.
#'
#' @return A numeric matrix with samples in rows and Roh_IS score in a column.
#'
#' @examples
#'
#' # Example does not matter as function will no be exported
compute_Roh_IS <- function(matches, RNA_tpm) {
  # Subset RNA_tpm
  sub_RNA_tpm <- RNA_tpm[matches, ]

  # Pseudocount of 0.01 for all genes
  sub_RNA_tpm <- sub_RNA_tpm + 0.01

  # Pseudocount of 1 for genes with 0 expr
  if (any(sub_RNA_tpm == 0)) sub_RNA_tpm[sub_RNA_tpm == 0] <- sub_RNA_tpm[sub_RNA_tpm == 0] + 1

  # Calculation: geometric mean (so-called log-average) [TPM, 0.01 offset]
  score <- apply(sub_RNA_tpm, 2, function(X) exp(mean(log(X))))

  return(data.frame(Roh_IS = score, check.names = FALSE))
}
