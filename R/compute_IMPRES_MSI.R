
#' Compute Immuno-Predictive Score (IMPRES) and Micro Satellite Instability (MSI) status score
#'
#' This function calculates IMPRES score by logical comparison of checkpoint gene pairs expression.
#'
#' @references
#'
#' Auslander,N.,Zhang,G.,Lee,J.S.,Frederick,D.T.,Miao,B.,Moll,T.,Tian,T.,Wei,Z., Madan, S.,
#' Sullivan, R.J., et al. (2018). Robust prediction of response to immune checkpoint blockade therapy in
#' metastatic melanoma. Nat. Med. 24, 1545–1549. https://doi.org/10.1038/s41591-018-0157-9.
#'
#' Fu, Y., Qi, L., Guo, W., Jin, L., Song, K., You, T., Zhang, S., Gu, Y., Zhao, W., and Guo, Z. (2019).
#' A qualitative transcriptional signature for predicting microsatellite instability status of right-sided
#' Colon Cancer. BMC Genomics 20, 769.
#'
#' @importFrom stats na.omit
#'
#' @param sig can be either 'IMPRES' or 'MSI'.
#' @param len the length of gene_1 vector.
#' @param match_F_1 numeric vector indicating the index of signature genes defined in 'gene_1' in `RNA_tpm`.
#' @param match_F_2 numeric vector indicating the index of signature genes defined in 'gene_2' in `RNA_tpm`.
#' @param RNA_tpm data.frame containing TPM values with HGNC symbols in rows and samples in columns.
#'
#' @return A numeric matrix with samples in rows and IMPRES score in a column.
#'
#' @examples
#'
#' # Example does not matter as function will no be exported
compute_IMPRES_MSI <- function(sig, len, match_F_1, match_F_2, RNA_tpm) {

  # Initialize variables
  F_pair_expr_A <- F_pair_expr_B <- SCORE_matrix <- matrix(0, len, ncol(RNA_tpm))
  colnames(SCORE_matrix) <- colnames(RNA_tpm)
  score <- vector("numeric", length = ncol(RNA_tpm))
  names(score) <- colnames(RNA_tpm)

  # Log2 transformation:
  log2_RNA_tpm <- as.data.frame(log2(RNA_tpm + 1))

  # Calculation:
  F_pair_expr_A <- log2_RNA_tpm[match_F_1, ]
  F_pair_expr_B <- log2_RNA_tpm[match_F_2, ]

  if (anyNA(F_pair_expr_A + F_pair_expr_B)) {
    remove_pairs <- as.vector(which(is.na(rowSums(F_pair_expr_A + F_pair_expr_B) == TRUE)))
  }

  SCORE_matrix <- F_pair_expr_A > F_pair_expr_B
  if (anyNA(SCORE_matrix)) {
    score <- colSums(SCORE_matrix, na.rm = TRUE)
    score <- (score * nrow(SCORE_matrix)) / (nrow(SCORE_matrix) - length(remove_pairs))
  } else {
    score <- colSums(SCORE_matrix)
  }

  df <- data.frame(score, check.names = FALSE)
  names(df)[1] <- sig

  return(df)
}
