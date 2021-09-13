#' Compute ligand-receptor pair weights from TPM bulk gene expression
#'
#' This function quantifies ligand-receptor interactions in the tumor microenvironment from
#' TPM bulk gene expression (Lapuente-Santana et al., bioRxiv, 2021), using prior knowledge
#' coming from ligand-receptor pair annotations from the database of (Ramilowski et al., Nat Commun, 2015).
#' Each ligand-receptor weight is defined as the minimum of the log2(TPM+1) expression of
#' the ligand and the receptor.
#'
#' @importFrom stats na.exclude
#' @importFrom utils head tail
#' @importFrom easierData get_intercell_networks get_group_lrpairs
#'
#' @param RNA_tpm A data.frame containing TPM values with HGNC symbols in rows and samples in columns.
#' @param cancer_type A string detailing the cancer type whose ligand-receptor pairs network will be used.
#' A pan-cancer network is selected by default, whose network represents the union of all
#' ligand-receptor pairs present across the 18 cancer types studied in (Lapuente-Santana et al., bioRxiv, 2021).
#' @param verbose A logical value indicating whether to display messages about the number of ligand-receptor
#' genes found in the gene expression data provided.
#'
#' @return A matrix of weights with samples in rows and ligand-receptor pairs in columns.
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
#' lrpair_weights[1:5, 1:5]
compute_LR_pairs <- function(RNA_tpm,
                             cancer_type = "pancan",
                             verbose = TRUE) {
  # Some checks
  if (is.null(RNA_tpm)) stop("TPM gene expression data not found")

  # Retrieve internal data
  intercell_networks <- suppressMessages(easierData::get_intercell_networks())
  group_lrpairs <- suppressMessages(easierData::get_group_lrpairs())

  # Gene expression data (log2 transformed)
  gene_expr <- log2(RNA_tpm + 1)
  genes <- rownames(gene_expr)

  # HGNC symbols are required
  if (any(grep("ENSG00000", genes))) stop("Hgnc gene symbols are required", call. = FALSE)

  gene_expr <- as.data.frame(gene_expr)

  # Cancer-specific LR pairs network
  intercell_network <- intercell_networks[[cancer_type]]
  LR_pairs <- unique(paste0(intercell_network$ligands, "_", intercell_network$receptors))

  # check what is the percentage of genes we have in our data
  all_lrpairs_genes <- unique(c(intercell_network$ligands, intercell_network$receptors))
  genes_kept <- intersect(rownames(gene_expr), all_lrpairs_genes)
  genes_left <- setdiff(all_lrpairs_genes, rownames(gene_expr))

  # check what is the percentage of regulated transcripts that we have in our data
  message("LR signature genes found in data set: ", length(genes_kept), "/", length(all_lrpairs_genes), " (", round(length(genes_kept) / length(all_lrpairs_genes), 3) * 100, "%)")

  # Compute L-R pairs
  LR_pairs_computed <- do.call(rbind, lapply(1:length(LR_pairs), function(x) {
    ligand <- sapply(strsplit(LR_pairs[x], split = "_", fixed = TRUE), head, 1)
    receptor <- sapply(strsplit(LR_pairs[x], split = "_", fixed = TRUE), tail, 1)

    pos_lr <- match(c(ligand, receptor), rownames(gene_expr))
    # When a ligand or receptor is not found, NA value should be returned.
    by_patient <- t(as.data.frame(apply(gene_expr[pos_lr, ], 2, min)))
    rownames(by_patient) <- LR_pairs[x]
    return(by_patient)
  }))
  LR_pairs_computed <- t(LR_pairs_computed)

  # Apply grouping to LRpairs data
  for (X in 1:length(group_lrpairs)) {
    keep <- unique(group_lrpairs[[X]]$main)
    remove <- unique(group_lrpairs[[X]]$involved_pairs)
    combo_name <- unique(group_lrpairs[[X]]$combo_name)

    pos_remove <- match(remove, colnames(LR_pairs_computed))
    pos_keep <- match(keep, colnames(LR_pairs_computed))

    colnames(LR_pairs_computed)[pos_keep] <- combo_name
    LR_pairs_computed <- LR_pairs_computed[, -pos_remove]
  }

  if (verbose) message("Ligand-Receptor pair weights computed \n")
  return(as.data.frame(LR_pairs_computed))
}
