#' Compute immune cell fractions from gene expression using quanTIseq
#'
#' Estimates cell fractions from TPM bulk gene expression
#' using quanTIseq method from Finotello et al., Genome Med, 2019.
#'
#' @references Finotello F, Mayer C, Plattner C, Laschober G, Rieder D,
#' Hackl H, Krogsdam A, Loncova Z, Posch W, Wilflingseder D, Sopper S,
#' Ijsselsteijn M, Brouwer TP, Johnson D, Xu Y, Wang Y, Sanders ME,
#' Estrada MV, Ericsson-Gonzalez P, Charoentong P, Balko J,
#' de Miranda NFDCC, Trajanoski Z. Molecular and pharmacological
#' modulators of the tumor immune contexture revealed by deconvolution
#' of RNA-seq data. Genome Medicine, 2019. 11(1):34.
#' https://doi.org/10.1186/s13073-019-0638-6
#'
#' @importFrom quantiseqr run_quantiseq
#'
#' @param RNA_tpm data.frame containing TPM values with HGNC symbols
#' in rows and samples in columns.
#' @param verbose logical value indicating whether to display messages
#' about the number of immune cell signature genes found in the gene
#' expression data provided.
#'
#' @return A numeric matrix of normalized enrichment scores
#' with samples in rows and cell types in columns.
#'
#' @export
#'
#' @examples
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
#' # Some genes are causing issues due to approved symbols matching more than one gene
#' genes_info <- reannotate_genes(cur_genes = rownames(RNA_tpm))
#' ## Remove non-approved symbols
#' non_na <- !is.na(genes_info$new_names)
#' RNA_tpm <- RNA_tpm[non_na, ]
#' genes_info <- genes_info[non_na, ]## Remove entries that are withdrawn
#' RNA_tpm <- RNA_tpm[-which(genes_info$new_names == "entry withdrawn"), ]
#' genes_info <- genes_info[-which(genes_info$new_names == "entry withdrawn"), ]
#' ## Identify duplicated new genes
#' newnames_dup <- unique(genes_info$new_names[duplicated(genes_info$new_names)])
#' newnames_dup_ind <- do.call(c, lapply(newnames_dup, function(X) which(genes_info$new_names == X)))
#' newnames_dup <- genes_info$new_names[newnames_dup_ind]
#' ## Retrieve data for duplicated genes
#' tmp <- RNA_tpm[genes_info$old_names[genes_info$new_names %in% newnames_dup],]
#' ## Remove data for duplicated genes
#' RNA_tpm <- RNA_tpm[-which(rownames(RNA_tpm) %in% rownames(tmp)),]
#' ## Aggregate data of duplicated genes
#' dup_genes <- genes_info$new_names[which(genes_info$new_names %in% newnames_dup)]
#' names(dup_genes) <- rownames(tmp)
#' if (anyDuplicated(newnames_dup)){
#' tmp2 <- stats::aggregate(tmp, by = list(dup_genes), FUN = "mean")
#' rownames(tmp2) <- tmp2$Group.1
#' tmp2$Group.1 <- NULL
#' }
#' # Put data together
#' RNA_tpm <- rbind(RNA_tpm, tmp2)
#'
#' # Computation of cell fractions (Finotello et al., Genome Med, 2019)
#' cell_fractions <- compute_cell_fractions(RNA_tpm = RNA_tpm)
compute_cell_fractions <- function(RNA_tpm = NULL,
                                   verbose = TRUE) {
  # Some checks
  if (is.null(RNA_tpm)) stop("TPM gene expression data not found")

  # HGNC symbols are required
  if (any(grep("ENSG00000", rownames(RNA_tpm)))) {
    stop("Hgnc gene symbols are required", call. = FALSE)
  }

  # Compute cell fractions: quanTIseq
  cell_fractions <- quantiseqr::run_quantiseq(
    expression_data = RNA_tpm, signature_matrix = "TIL10",
    is_arraydata = FALSE, is_tumordata = TRUE, scale_mRNA = TRUE
  )
  cell_fractions$Sample <- NULL

  # Change names to match opt models
  new_cellnames <- c(
    "B", "M1", "M2", "Monocyte", "Neutrophil",
    "NK", "CD4 T", "CD8+ T", "Treg", "DC", "Other"
  )
  colnames(cell_fractions) <- new_cellnames

  # Fix estimation issue with Tregs and CD4 T
  cell_fractions[, "CD4 T"] <- cell_fractions[, "CD4 T"] + cell_fractions[, "Treg"]

  if (verbose) message("Cell fractions computed! \n")
  return(cell_fractions)
}
