#' Compute cell-cell pair score
#'
#' This function derives a score for each cell-cell pair feature.
#'
#' @param celltype1 string character with first cell type involved in the interaction.
#' @param celltype2 string character with second cell type involved in the interaction.
#' @param intercell_network matrix with data on cell types interaction network.
#' @param lrpairs_binary binary vector displaying LR pairs with non-zero frequency.
#' @param lr_frequency numeric vector with LR pairs frequency across the whole TCGA database.
#' @param compute_log boolean variable in order to take the log of the weighted score.
#'
#' @return A numeric vector with weighted scores.
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
#' # remove ligand receptor pairs that are always NA
#' na_lrpairs <- apply(lrpair_weights, 2, function(x) {
#'   all(is.na(x))
#' })
#' lrpair_weights <- lrpair_weights[, na_lrpairs == FALSE]
#'
#' # binarize the data: set a threshold to 10 TPM,
#' # only pairs where both ligand and receptor have TPM > 10 are kept
#' lrpairs_binary <- ifelse(lrpair_weights > log2(10 + 1), 1, 0)
#'
#' # keep only the LR.pairs for which I have (non-zero) frequencies in the TCGA
#' lr_frequency <- easierdata_eh[["EH6684"]]
#' lrpairs_binary <- lrpairs_binary[, colnames(lrpairs_binary) %in% names(lr_frequency)]
#'
#' # cancer type specific network
#' intercell_networks <- easierdata_eh[["EH6683"]]
#' intercell_network <- intercell_networks[["pancan"]]
#' celltypes <- unique(c(as.character(intercell_network$cell1), as.character(intercell_network$cell2)))
#' celltype1 <- celltypes[1]
#' celltype2 <- celltypes[1]
#'
#' # compute the CC score for each patient
#' CCpair_score <- compute_CCpair_score(celltype1, celltype2, intercell_network,
#'   lrpairs_binary, lr_frequency,
#'   compute_log = TRUE
#' )
compute_CCpair_score <- function(celltype1,
                                 celltype2,
                                 intercell_network,
                                 lrpairs_binary,
                                 lr_frequency,
                                 compute_log = TRUE) {

  # consider the LR interactions between the two cell types
  CC_network <- intercell_network[intersect(which(intercell_network$cell1 == celltype1), which(intercell_network$cell2 == celltype2)), ]
  CC_LRpairs <- paste(CC_network$ligands, CC_network$receptors, sep = "_")

  # extract the corresponding data for all patients
  ix <- match(CC_LRpairs, colnames(lrpairs_binary))
  CC_LR_data <- lrpairs_binary[, ix[!is.na(ix)]]

  # and the LR frequecies
  CC_lr_frequency <- lr_frequency[colnames(CC_LR_data)]

  # multiply each row of the matrix (i.e. each patient data) for the vector with the frequencies
  CC_LR_data_weighted <- t(t(CC_LR_data) * 1 / CC_lr_frequency)

  # compute the cell cell interaction score as the sum of the LR weighted pairs
  CC_score <- apply(CC_LR_data_weighted, 1, sum)

  # if we use the weighted score taking the log might be better
  if (compute_log == TRUE) {
    CC_score <- log2(CC_score + 1)
  }
}
