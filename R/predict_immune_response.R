#' Compute predicted immune response
#'
#' This function calculates predictions of patients' immune response using the quantitative descriptors
#' data as input features and the model parameters defined during training.
#'
#' @importFrom utils combn
#' @importFrom stats na.omit
#' @importFrom BiocParallel register bplapply MulticoreParam
#' @importFrom easierData get_opt_models get_opt_xtrain_stats
#'
#' @export
#'
#' @param pathways numeric matrix with pathways activity (rows = samples; columns = pathways).
#' @param immunecells numeric matrix with immune cell quantification (rows = samples; columns = cell types).
#' @param tfs numeric matrix with transcription factors activity (rows = samples; columns = transcription factors).
#' @param lrpairs numeric matrix with ligand-receptor weights (rows = samples; columns = ligand-receptor pairs).
#' @param ccpairs numeric matrix with cell-cell scores (rows = samples; columns = cell-cell pairs).
#' @param cancer_type character string indicating which cancer-specific model should be used to compute the predictions.
#' @param verbose logical flag indicating whether to display messages about the process.
#'
#' @return A list containing the predictions for each quantitative descriptor and for each task.
#' Given that the model training was repeated 100 times with randomized-cross validation, a set of 100 predictions is returned.
#'
#' @examples
#' # Load exemplary dataset (Mariathasan et al., Nature, 2018) from easierData.
#' # Original processed data is available from IMvigor210CoreBiologies package.
#' library("easierData")
#' dataset_mariathasan <- easierData::get_Mariathasan2018_PDL1_treatment()
#' RNA_tpm <- dataset_mariathasan@assays@data@listData[["tpm"]]
#' RNA_counts <- dataset_mariathasan@assays@data@listData[["counts"]]
#' cancer_type <- dataset_mariathasan@metadata$cancertype
#'
#' # Computation of cell fractions (Finotello et al., Genome Med, 2019)
#' cell_fractions <- compute_cell_fractions(RNA_tpm = RNA_tpm)
#'
#' # Computation of pathway scores (Holland et al., BBAGRM, 2019; Schubert et al., Nat Commun, 2018)
#' pathway_activity <- compute_pathway_activity(
#'   RNA_counts = RNA_counts,
#'   remove_sig_genes_immune_response = TRUE
#' )
#'
#' # Computation of TF activity (Garcia-Alonso et al., Genome Res, 2019)
#' tf_activity <- compute_TF_activity(
#'   RNA_tpm = RNA_tpm
#' )
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
#'
#' # Predict patients' immune response
#' predictions_immune_response <- predict_immune_response(
#'   pathways = pathway_activity,
#'   immunecells = cell_fractions,
#'   tfs = tf_activity,
#'   lrpairs = lrpair_weights,
#'   ccpairs = ccpair_scores,
#'   cancer_type = cancer_type
#' )
predict_immune_response <- function(pathways = NULL,
                                    immunecells = NULL,
                                    tfs = NULL,
                                    lrpairs = NULL,
                                    ccpairs = NULL,
                                    cancer_type,
                                    verbose = TRUE) {
  if (missing(cancer_type)) stop("cancer type needs to be specified")
  if (all(is.null(pathways), is.null(immunecells), is.null(tfs), is.null(lrpairs), is.null(ccpairs))) stop("none signature specified")

  # Initialize variables
  views <- c(
    pathways = "gaussian",
    immunecells = "gaussian",
    tfs = "gaussian",
    lrpairs = "gaussian",
    ccpairs = "gaussian"
  )

  # Check which views are missing
  miss_views <- c(
    ifelse(missing(pathways), NA, 1),
    ifelse(missing(immunecells), NA, 2),
    ifelse(missing(tfs), NA, 3),
    ifelse(missing(lrpairs), NA, 4),
    ifelse(missing(ccpairs), NA, 5)
  )

  # Single views
  view_simples <- lapply(miss_views[!is.na(miss_views)], function(X) {
    tmp <- views[X]
    return(tmp)
  })

  # All corresponding views
  view_combinations <- view_simples

  # Retrieve internal data
  opt_models <- suppressMessages(easierData::get_opt_models())
  opt_xtrain_stats <- suppressMessages(easierData::get_opt_xtrain_stats())

  compute_prediction <- function(view, verbose, opt_models, opt_xtrain_stats, cancer_type) {
    view_info <- view_combinations[[view]]
    view_name <- paste(names(view_info), collapse = "_")
    view_data <- lapply(tolower(names(view_info)), function(x) as.data.frame(get(x)))
    names(view_data) <- names(view_info)
    if (verbose) message("Computing predictions using ", view_name, "...\n")

    # Initialize variables
    opt_model_cancer_view_spec <- lapply(view_name, function(X) {
      return(opt_models[[cancer_type]][[X]])
    })
    names(opt_model_cancer_view_spec) <- view_name
    opt_xtrain_stats_cancer_view_spec <- lapply(view_name, function(X) {
      return(opt_xtrain_stats[[cancer_type]][[X]])
    })
    names(opt_xtrain_stats_cancer_view_spec) <- view_name

    # Predict immune response using RMTLR model parameters
    prediction_view <- predict_with_rmtlr(
      view_name = view_name,
      view_info = view_info,
      view_data = view_data,
      opt_model_cancer_view_spec = opt_model_cancer_view_spec,
      opt_xtrain_stats_cancer_view_spec = opt_xtrain_stats_cancer_view_spec
    )

    return(prediction_view)
  }
  # Parallelize views model predictions
  BiocParallel::register(BiocParallel::MulticoreParam(workers = 2))
  all_predictions <- BiocParallel::bplapply(1:length(view_combinations), FUN = compute_prediction, verbose = verbose, opt_models, opt_xtrain_stats, cancer_type)

  names(all_predictions) <- sapply(1:length(view_combinations), function(X) {
    paste(names(view_combinations[[X]]), collapse = "_")
  })
  return(all_predictions)
}
