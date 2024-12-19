#' Compute random scores of the immune resistance
#' program used in the computation of repressed
#' immune resistance signature (RIR) score.
#'
#' Calculates random scores to yield a robust estimate of
#' the immune resistance program values. This is used
#' by get_OE_bulk function.
#'
#' The source code was provided by original work:
#' https://github.com/livnatje/ImmuneResistance
#'
#' @references Jerby-Arnon, L., Shah, P., Cuoco, M.S.,
#' Rodman, C., Su, M.-J., Melms, J.C., Leeson, R., Kanodia,
#' A., Mei, S., Lin, J.-R., et al. (2018). A Cancer Cell
#' Program Promotes T Cell Exclusion and Resistance to
#' Checkpoint Blockade. Cell 175, 984–997.e24.
#' https://doi.org/10.1016/j.cell.2018.09.006.
#'
#' @param r list containing a numeric matrix with bulk RNA-Seq
#' data (tpm values) and a character string with the available
#' gene names.
#' @param genes_dist_q factor variable obtained as output from
#' the function discretize. Original work binned genes into 50
#' expression bins according their average gene expression
#' across samples.
#' @param b_sign logical vector representing whether signature
#' genes were found in bulk tpm matrix.
#' @param num_rounds integer value related to the number of random
#' gene signatures samples to be computed for normalization. Original
#' work indicates that 1000 random signatures were sufficient to yield
#' an estimate of the expected value.
#' @param full_flag logical flag indicating whether to return also
#' random scores.
#' @param random_seed integer value to set a seed for the selection
#' of random genes used to generate a random score.
#'
#' @return A numeric vector containing the estimated random score
#' for each sample.
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
#' # Log2 transformation:
#' log2_RNA_tpm <- log2(RNA_tpm + 1)
#'
#' # Prepare input data
#' r <- list()
#' r$tpm <- log2_RNA_tpm
#' r$genes <- rownames(log2_RNA_tpm)
#'
#' # Gene signature of immune resistance program
#' score_signature_genes <- suppressMessages(easierData::get_scores_signature_genes())
#' RIR_gene_signature <- score_signature_genes$RIR
#'
#' # Compute gene average expression across samples
#' r$genes_dist <- r$genes_mean <- rowMeans(r$tpm)
#'
#' # Center gene expression matrix
#' r$zscores <- sweep(r$tpm, 1, r$genes_mean, FUN = "-")
#'
#' # Bin genes into 50 expression bins according to their average
#' r$genes_dist_q <- easier:::discretize(r$genes_dist, n_cat = 50)
#'
#' # Match genes from exc.down signature with genes from expression matrix
#' b_sign <- is.element(r$genes, RIR_gene_signature[["exc.down"]])
#'
#' # Compute random score:
#' rand_scores <- easier:::get_semi_random_OE(r, r$genes_dist_q, b_sign)
#' }
get_semi_random_OE <- function(r,
                               genes_dist_q,
                               b_sign,
                               num_rounds = 1000,
                               full_flag = FALSE,
                               random_seed = 1234) {
  # Previous name: get.random.sig.scores
  sign_q <- as.matrix(table(genes_dist_q[b_sign]))
  q <- rownames(sign_q)
  idx_all <- c()
  B <- matrix(
    data = FALSE, nrow = length(genes_dist_q),
    ncol = num_rounds
  )
  Q <- matrix(
    data = 0, nrow = length(genes_dist_q),
    ncol = num_rounds
  )
  for (i in seq_len(nrow(sign_q))) {
    num_genes <- sign_q[i]
    if (num_genes > 0) {
      idx <- which(is.element(genes_dist_q, q[i]))
      for (j in seq_len(num_rounds)) {
        idxj <- sample(idx, num_genes)
        Q[i, j] <- sum(B[idxj, j] == TRUE)
        B[idxj, j] <- TRUE
      }
    }
  }
  rand_scores <- apply(B, 2, function(x) colMeans(r$zscores[x, ]))
  if (full_flag) {
    return(rand_scores)
  }
  rand_scores <- rowMeans(rand_scores)
  return(rand_scores)
}
