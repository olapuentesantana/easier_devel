#' Compute repressed immune resistance signature (RIR) score
#'
#' This function calculates RIR score defined by combining a set of gene signatures associated
#' with downregulation of T cell exclusion, post-treatment and functional resistance. These gene
#' signatures were obtained from:
#' https://github.com/livnatje/ImmuneResistance/blob/master/Results/Signatures/resistance.program.RData.
#'
#' @references Jerby-Arnon, L., Shah, P., Cuoco, M.S., Rodman, C., Su, M.-J., Melms, J.C., Leeson, R., Kanodia, A., Mei, S., Lin, J.-R., et al. (2018).
#' A Cancer Cell Program Promotes T Cell Exclusion and Resistance to Checkpoint Blockade. Cell 175, 984–997.e24. https://doi.org/10.1016/j.cell.2018.09.006.
#'
#' @importFrom stats na.omit
#'
#' @param RNA_tpm data.frame containing TPM values with HGNC symbols in rows and samples in columns.
#' @param RIR_program list with gene signatures included in the immune resistance program from Jerby-Arnon et al., 2018.
#'
#' @return A numeric matrix with samples in rows and RIR score in a column.
#'
#' @examples
#'
#' # Example does not matter as function will no be exported
compute_RIR <- function(RNA_tpm,
                        RIR_program) {

  # Log2 transformation:
  log2_RNA_tpm <- log2(RNA_tpm + 1)

  # Prepare input data
  r <- list()
  r$tpm <- log2_RNA_tpm
  r$genes <- rownames(log2_RNA_tpm)

  # Apply function to calculate OE:
  res_scores <- get_OE_bulk(r, gene_sign = RIR_program, verbose = TRUE)

  # Merge as recommend by authors
  res <- cbind.data.frame(
    excF.up = rowMeans(res_scores[, c("exc.up", "exc.seed.up")]),
    excF.down = rowMeans(res_scores[, c("exc.down", "exc.seed.down")]),
    res.up = rowMeans(res_scores[, c("trt.up", "exc.up", "exc.seed.up")]),
    res.down = rowMeans(res_scores[, c("trt.down", "exc.down", "exc.seed.down")]),
    res_scores
  )

  res <- cbind.data.frame(
    resF.up = res[, "res.up"] + res[, "fnc.up"],
    resF.down = res[, "res.down"] + res[, "fnc.down"],
    res
  )

  # Keep that signature considered to be relevant
  keep_sig <- c("resF.down")
  score <- as.matrix(res[, colnames(res) %in% keep_sig])
  rownames(score) <- colnames(log2_RNA_tpm)

  return(data.frame(RIR = score, check.names = FALSE))
}
