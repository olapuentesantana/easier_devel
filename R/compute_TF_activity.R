#' Compute transcription factor activity from gene expression using DoRothEA
#'
#' This function infers transcription factor activity from TPM bulk gene expression
#' using DoRothEA method (Garcia-Alonso et al., Genome Res, 2019).
#'
#' @importFrom dorothea run_viper
#' @importFrom stats na.exclude
#' @importFrom dplyr filter
#' @importFrom easierData get_TCGA_mean_pancancer get_TCGA_sd_pancancer
#'
#' @param RNA_tpm data.frame containing TPM values with HGNC symbols in rows and samples in columns.
#' @param verbose logical value indicating whether to display messages about the number of regulated
#' genes found in the gene expression data provided.
#'
#' @return A numeric matrix of activity scores with samples in rows and transcription factors in columns.
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
#' # Computation of TF activity (Garcia-Alonso et al., Genome Res, 2019)
#' tf_activity <- compute_TF_activity(
#'   RNA_tpm = RNA_tpm
#' )
compute_TF_activity <- function(RNA_tpm,
                                verbose = TRUE) {
  # Some checks
  if (is.null(RNA_tpm)) stop("TPM gene expression data not found")

  # Retrieve internal data
  TCGA_mean_pancancer <- suppressMessages(easierData::get_TCGA_mean_pancancer())
  TCGA_sd_pancancer <- suppressMessages(easierData::get_TCGA_sd_pancancer())

  # Gene expression data
  tpm <- RNA_tpm
  genes <- rownames(tpm)

  # HGNC symbols are required
  if (any(grep("ENSG00000", genes))) stop("hgnc gene symbols are required", call. = FALSE)

  # Log transformed expression matrix (log2[tpm+1]): expression matrix scaled and recentered.
  gene_expr <- calc_z_score(t(tpm), mean = TCGA_mean_pancancer, sd = TCGA_sd_pancancer)

  # redefine gene names to match transcripts for viper
  E <- t(gene_expr)
  newNames <- sapply(rownames(E), function(x) {
    # strsplit(x, "\\.")[[1]][1]
    zz_tmp <- strsplit(x, "\\.")[[1]]
    paste0(zz_tmp[1:(length(zz_tmp) - 1)], collapse = "-")
  })
  rownames(E) <- newNames

  # data extracted from publication
  regulons <- dplyr::filter(dorothea::dorothea_hs, .data$confidence %in% c("A", "B"))
  all_regulated_transcripts <- unique(regulons$target)
  all_tfs <- unique(regulons$tf)

  # check what is the percentage of genes we have in our data
  genes_kept <- intersect(rownames(E), all_regulated_transcripts)
  genes_left <- setdiff(all_regulated_transcripts, rownames(E))

  # check what is the percentage of regulated transcripts that we have in our data
  if (verbose) {
    message(
      "Regulated transcripts found in data set: ", length(genes_kept), "/",
      length(all_regulated_transcripts), " (",
      round(length(genes_kept) / length(all_regulated_transcripts), 3) * 100, "%)"
    )
  }

  # TF activity: run viper
  tf_activity <- dorothea::run_viper(
    input = E, regulons = regulons,
    options = list(method = "none", minsize = 4, eset.filter = FALSE, cores = 1, verbose = FALSE)
  )

  # Samples as rows, TFs as columns
  tf_activity <- t(tf_activity)

  if (verbose) message("TF activity computed! \n")

  return(as.data.frame(tf_activity))
}
