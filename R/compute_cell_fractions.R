#' Compute immune cell fractions from gene expression using quanTIseq
#'
#' This function estimates cell fractions from TPM bulk gene expression
#' using quanTIseq method from (Finotello et al., Genome Med, 2019).
#'
#' @importFrom quantiseqr run_quantiseq
#'
#' @param RNA_tpm data.frame containing TPM values with HGNC symbols in rows and samples in columns.
#' @param verbose logical value indicating whether to display messages about the number of immune cell
#' signature genes found in the gene expression data provided.
#'
#' @return A numeric matrix of normalized enrichment scores with samples in rows and cell types in columns.
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
#' # Computation of cell fractions (Finotello et al., Genome Med, 2019)
#' cell_fractions <- compute_cell_fractions(RNA_tpm = RNA_tpm)
compute_cell_fractions <- function(RNA_tpm,
                                   verbose = TRUE
) {
  # Some checks
  if (is.null(RNA_tpm)) stop("TPM gene expression data not found")

  # HGNC symbols are required
  if (any(grep("ENSG00000", rownames(RNA_tpm)))) stop("Hgnc gene symbols are required", call. = FALSE)

  # Cell fractions: run deconvolute
  cell_fractions <- quantiseqr::run_quantiseq(
    expression_data = RNA_tpm, signature_matrix = "TIL10",
    is_arraydata = FALSE, is_tumordata = TRUE, scale_mRNA = TRUE
  )

  cell_fractions$Sample <- NULL
  # Samples as rows, immune cells as columns
  old_cellnames <- colnames(cell_fractions)
  new_cellnames <- c("B", "M1", "M2", "Monocyte", "Neutrophil", "NK", "CD4 T", "CD8+ T", "Treg", "DC", "Other")
  colnames(cell_fractions) <- new_cellnames
  cell_fractions[, "CD4 T"] <- cell_fractions[, "CD4 T"] + cell_fractions[, "Treg"]

  if (verbose) message("Cell fractions computed! \n")
  return(cell_fractions)
}
