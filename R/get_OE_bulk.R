#' Compute overall expression (OE) of the immune
#' resistance program used in the computation of
#' repressed immune resistance signature (RIR) score.
#'
#' This function calculates the overall expression of
#' the immune resistance program which is based on a set
#' of gene signatures associated with T cell exclusion,
#' post-treatment and functional resistance.
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
#' @param gene_sign list containing different character strings
#' associated with subsets of the resistance program.
#' @param num_rounds integer value related to the number of random
#' gene signatures samples to be computed for normalization.
#' Original work indicates that 1000 random signatures were sufficient
#' to yield an estimate of the expected value.
#' @param full_flag logical flag indicating whether to return also
#' random scores.
#' @param verbose logical flag indicating whether to display messages
#' about the process.
#'
#' @return A numeric matrix with computed scores for each sample and
#' subset of signatures included in the immune resistance program
#' (rows = samples; columns = gene signatures)
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
#' # Apply function to calculate OE:
#' res_scores <- easier:::get_OE_bulk(r, gene_sign = RIR_gene_signature, verbose = TRUE)
#' }
get_OE_bulk <- function(r,
                        gene_sign = NULL,
                        num_rounds = 1000,
                        full_flag = FALSE,
                        verbose = TRUE) {
  r$genes_mean <- rowMeans(r$tpm)
  r$zscores <- sweep(r$tpm, 1, r$genes_mean, FUN = "-")
  r$genes_dist <- r$genes_mean
  r$genes_dist_q <- discretize(r$genes_dist, n_cat = 50)

  r$sig_scores <- matrix(data = 0, nrow = ncol(r$tpm), ncol = length(gene_sign))
  sig_names <- names(gene_sign)
  colnames(r$sig_scores) <- sig_names
  r$sig_scores_raw <- r$sig_scores
  rand_flag <- is.null(r$rand_scores) | !all(is.element(
    names(gene_sign),
    colnames(r$rand_scores)
  ))
  if (rand_flag) {
    # if (verbose) message("Computing also random scores...", "\n")
    r$rand_scores <- r$sig_scores
  }
  for (i in sig_names) {
    b_sign <- is.element(r$genes, gene_sign[[i]])
    if (sum(b_sign) < 2) {
      next()
    }
    if (rand_flag) {
      rand_scores <- get_semi_random_OE(r, r$genes_dist_q, b_sign,
        num_rounds = num_rounds, random_seed = 1234
      )
    } else {
      rand_scores <- r$rand_scores[, i]
    }
    raw_scores <- colMeans(r$zscores[b_sign, ])
    final_scores <- raw_scores - rand_scores
    r$sig_scores[, i] <- final_scores
    r$sig_scores_raw[, i] <- raw_scores
    r$rand_scores[, i] <- rand_scores
  }
  if (full_flag) {
    return(r)
  }
  sig_scores <- r$sig_scores
  return(sig_scores)
}
