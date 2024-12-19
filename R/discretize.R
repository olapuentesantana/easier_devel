#' Converts a continuous variable into categorical
#'
#' It is used to bin continuous gene expression values
#' from a given gene signature into categories.
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
#' @importFrom stats quantile
#'
#' @param v numeric vector with gene mean expression across samples.
#' @param n_cat number of categories to bin continuous values,
#' here gene expression values.
#'
#' @return A numeric vector providing an integer value (e.g. category)
#' for each gene.
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
#'
#' # Gene signature of immune resistance program
#' score_signature_genes <- suppressMessages(easierData::get_scores_signature_genes())
#' RIR_gene_signature <- score_signature_genes$RIR
#'
#' # Compute gene average expression across samples
#' r$genes_dist <- r$genes_mean <- rowMeans(r$tpm)
#'
#' # Bin genes into 50 expression bins according to their average
#' r$genes_dist_q <- easier:::discretize(r$genes_dist, n_cat = 50)
#' }
discretize <- function(v, n_cat) {
  q1 <- stats::quantile(v, seq(from = (1 / n_cat), to = 1, by = (1 / n_cat)))
  u <- matrix(nrow = length(v))
  for (i in 2:n_cat) {
    u[(v >= q1[i - 1]) & (v < q1[i])] <- i
  }
  return(u)
}
