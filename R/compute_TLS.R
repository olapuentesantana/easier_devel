#' Computation of tertiary lymphoid structures signature (TLS) score
#'
#' This function calculates TLS score as the geometric-mean of the expression of its signature genes.
#'
#' @references Cabrita, R., Lauss, M., Sanna, A., Donia, M., Skaarup Larsen, M., Mitra, S., Johansson, I., Phung, B.,
#' Harbst, K., Vallon-Christersson, J., et al. (2020). Tertiary lymphoid structures improve immunotherapy and survival
#' in melanoma. Nature 577, 561–565.
#'
#' @importFrom stats na.omit
#'
#' @param matches numeric vector indicating the index of signature genes in `RNA_tpm`.
#' @param RNA_tpm data.frame containing TPM values with HGNC symbols in rows and samples in columns.
#'
#' @return A numeric matrix with samples in rows and TLS score in a column.
#'
#' @examples
#'
#' # Example does not matter as function will no be exported
compute_TLS <- function(matches, RNA_tpm) {
  # Subset RNA_tpm
  sub_RNA_tpm <- RNA_tpm[matches, ]

  # Calculation: geometric mean (so-called log-average) [TPM, 1 offset]
  geom_mean <- apply(sub_RNA_tpm, 2, function(X) exp(mean(log2(X + 1))))

  return(data.frame(TLS = geom_mean, check.names = FALSE))
}
