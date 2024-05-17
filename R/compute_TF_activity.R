#' Compute transcription factor activity from gene
#' expression using DoRothEA
#'
#' Infers transcription factor (TF) activity from TPM bulk
#' gene expression using DoRothEA method from Garcia-Alonso et al.,
#' Genome Res, 2019.
#'
#' @references Garcia-Alonso L, Holland CH, Ibrahim MM, Turei D,
#' Saez-Rodriguez J. "Benchmark and integration of resources for
#' the estimation of human transcription factor activities."
#' Genome Research. 2019. DOI: 10.1101/gr.240663.118.
#'
#' @import magrittr
#' @importFrom stats na.exclude
#' @importFrom dplyr filter
#' @importFrom decoupleR get_collectri get_dorothea run_viper
#' @importFrom tidyr pivot_wider
#' @importFrom tibble column_to_rownames
#' @importFrom easierData get_TCGA_mean_pancancer get_TCGA_sd_pancancer
#'
#' @param RNA_tpm data.frame containing TPM values with HGNC symbols
#' in rows and samples in columns.
#' @param regulon_net string indicating the regulon network to be used.
#' @param verbose logical value indicating whether to display messages
#' about the number of regulated
#' genes found in the gene expression data provided.
#'
#' @return A numeric matrix of activity scores with samples in rows
#' and TFs in columns.
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
#' # Computation of TF activity (Garcia-Alonso et al., Genome Res, 2019)
#' tf_activity <- compute_TF_activity(
#'   RNA_tpm = RNA_tpm,
#'   regulon_net = "collectri"
#' )
compute_TF_activity <- function(RNA_tpm = NULL,
                                regulon_net = c("collectri", "dorothea"),
                                verbose = TRUE) {
  # Some checks
  if (is.null(RNA_tpm)) stop("TPM gene expression data not found")
  if (length(regulon_net) > 1) stop("Please select just one regulon network")
  
  # Retrieve internal data
  TCGA_mean_pancancer <- suppressMessages(easierData::get_TCGA_mean_pancancer())
  TCGA_sd_pancancer <- suppressMessages(easierData::get_TCGA_sd_pancancer())

  # Gene expression data
  tpm <- RNA_tpm
  genes <- rownames(tpm)

  # HGNC symbols are required
  if (any(grep("ENSG00000", genes))) stop("hgnc gene symbols are required", call. = FALSE)
  
  # Log transformed expression matrix (log2[tpm+1]):
  # expression matrix scaled and recentered.
  gene_expr <- calc_z_score(t(tpm),
                            mean = TCGA_mean_pancancer,
                            sd = TCGA_sd_pancancer
  )

  # redefine gene names to match TF-target network
  E <- t(gene_expr)
  newNames <- gsub(".", "-", rownames(E), fixed = TRUE)
  rownames(E) <- newNames

  # collectTRI network
  if(regulon_net == "collectri"){
    
    net <- decoupleR::get_collectri(organism='human', split_complexes=FALSE)
  
  # dorothea network
  }else if(regulon_net == "dorothea"){
    
    net <- decoupleR::get_dorothea(organism='human', 
                                   levels = c("A", "B"))
    
  }

  all_regulated_transcripts <- unique(net$target)

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

  # Remove genes with all NA/Inf values
  E <- E[!is.na(apply(E, 1, sum)), ]
  E <- E[!is.infinite(apply(E, 1, sum)), ]

  # TF activity: run viper
  tf_activity_df <- decoupleR::run_viper(mat = E,
                                         network = net,
                                         .source='source',
                                         .target='target',
                                         .mor='mor',
                                         minsize = 4,
                                         eset.filter = FALSE,
                                         method = "none")
  # To matrix
  tf_activity <- tf_activity_df %>%
    tidyr::pivot_wider(id_cols = condition,
                       names_from = source,
                       values_from = score) %>%
    tibble::column_to_rownames("condition")

  if (verbose) message("TF activity computed! \n")

  return(as.data.frame(tf_activity))
}
