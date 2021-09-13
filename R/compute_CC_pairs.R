#' Compute cell-cell interactions scores using computed ligand-receptor weights
#'
#' This function scores cell-cell interactions in the tumor microenvironment using
#' ligand-receptor weights as input (Lapuente-Santana et al., bioRxiv, 2021).
#'
#' @importFrom easierData get_lr_frequency_TCGA get_intercell_networks
#'
#' @param lrpairs output of the compute_LR_pairs function. A matrix of log2(TPM +1) weights with
#' with samples in rows and ligand-receptor pairs in columns.
#' @param cancer_type string detailing the cancer type whose cell-cell interaction network will be used.
#' A pan-cancer network is selected by default, whose network represents the union of all
#' ligand-receptor pairs present across the 18 cancer types studied in (Lapuente-Santana et al., bioRxiv, 2021).
#' @param verbose logical value indicating whether to display informative messages about the process.
#'
#' @return A matrix of scores with samples in rows and cell-cell pairs in columns.
#'
#' @export
#'
#' @examples
#' # Load exemplary dataset (Mariathasan et al., Nature, 2018) from easierData.
#' # Original processed data is available from IMvigor210CoreBiologies package.
#' library("easierData")
#' dataset_mariathasan <- easierData::get_Mariathasan2018_PDL1_treatment()
#' RNA_tpm <- dataset_mariathasan@assays@data@listData[["tpm"]]
#'
#' # Computation of ligand-receptor pair weights
#' lrpair_weights <- compute_LR_pairs(
#'   RNA_tpm = RNA_tpm,
#'   cancer_type = "pancan"
#' )
#'
#' # Computation of cell-cell interaction scores
#' ccpair_scores <- compute_CC_pairs(
#'   lrpairs = lrpair_weights,
#'   cancer_type = "pancan"
#' )
compute_CC_pairs <- function(lrpairs,
                             cancer_type = "pancan",
                             verbose = TRUE) {
  # Some checks
  if (is.null(lrpairs)) stop("ligand-receptor pair weights not found")

  # Retrieve internal data
  lr_frequency <- suppressMessages(easierData::get_lr_frequency_TCGA())
  intercell_networks <- suppressMessages(easierData::get_intercell_networks())

  # remove ligand receptor pairs that are always NA
  na_lrpairs <- apply(lrpairs, 2, function(x) {
    all(is.na(x))
  })
  lrpairs <- lrpairs[, na_lrpairs == FALSE]

  # binarize the data: set a threshold to 10 TPM, only pairs where both ligand and receptor have
  # TPM > 10 are kept
  lrpairs_binary <- ifelse(lrpairs > log2(10 + 1), 1, 0)
  # sum(LR.pairs.binary)/(ncol(LR.pairs.binary)*nrow(LR.pairs.binary))*100 #percentage of

  # keep only the LR.pairs for which I have (non-zero) frequencies in the TCGA
  lrpairs_binary <- lrpairs_binary[, colnames(lrpairs_binary) %in% names(lr_frequency)]

  # cancer type specific network
  intercell_network <- intercell_networks[[cancer_type]]

  # compute the CC score for each patient
  celltypes <- unique(c(as.character(intercell_network$cell1), as.character(intercell_network$cell2)))

  CC_pairs_score <- do.call(cbind, lapply(celltypes, function(celltype1) {
    do.call(cbind, lapply(celltypes, function(celltype2) {
      compute_CCpair_score(celltype1, celltype2, intercell_network,
        lrpairs_binary, lr_frequency,
        compute_log = TRUE
      )
    }))
  }))

  metadata_CC_pairs <- do.call(rbind, lapply(celltypes, function(celltype1) {
    do.call(rbind, lapply(celltypes, function(celltype2) {
      data.frame(
        CCpair = gsub(" ", "", paste(celltype1, celltype2, sep = "_")),
        celltype1 = celltype1, celltype2 = celltype2
      )
    }))
  }))

  colnames(CC_pairs_score) <- metadata_CC_pairs$CCpair

  if (verbose) message("CC pairs computed \n")
  return(as.data.frame(CC_pairs_score))
}
