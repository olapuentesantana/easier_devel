#' Collects Regularized Multi-Task Linear Regression (RMTLR) cancer-specfic model predictions
#'
#' This function puts together predicted immune response by using a cancer-specifc model learned
#' from training with RMTLR algorithm.
#'
#' @importFrom stats na.omit
#'
#' @export
#'
#' @param view_name character string containing the name of the input view.
#' @param view_info character string informing about the family of the input data.
#' @param view_data list containing the data for each input view.
#' @param opt_model_cancer_view_spec cancer-view-specific model feature parameters learned during training.
#' @param opt_xtrain_stats_cancer_view_spec cancer-view-specific features mean and standard deviation of the training set.
#' @param verbose logical flag indicating whether to display messages about the process.
#'
#' @return A lList of predictions matrices, one for each tasks (rows = samples; columns = [runs).
#'
#' @examples
#' # Load exemplary dataset (Mariathasan et al., Nature, 2018) from easierData.
#' # Original processed data is available from IMvigor210CoreBiologies package.
#' library("easierData")
#' dataset_mariathasan <- easierData::get_Mariathasan2018_PDL1_treatment()
#' RNA_tpm <- dataset_mariathasan@assays@data@listData[["tpm"]]
#' cancer_type <- dataset_mariathasan@metadata$cancertype
#'
#' # Computation of cell fractions
#' cell_fractions <- compute_cell_fractions(RNA_tpm = RNA_tpm)
#'
#' view_name <- "immunecells"
#' view_info <- c(immunecells = "gaussian")
#' view_data <- list(immunecells = as.data.frame(cell_fractions))
#'
#' # Predict using rmtlr
#' prediction_view <- predict_with_rmtlr(
#'   view_name = view_name,
#'   view_info = view_info,
#'   view_data = view_data,
#'   cancer_type = cancer_type
#' )
predict_with_rmtlr <- function(view_name,
                               view_info,
                               view_data,
                               opt_model_cancer_view_spec,
                               opt_xtrain_stats_cancer_view_spec,
                               verbose = TRUE) {

  # initialize variables
  P <- length(view_info)
  K <- ncol(opt_model_cancer_view_spec[[1]][[1]])
  standardize_any <- TRUE
  tasks <- names(opt_model_cancer_view_spec[[1]])

  # Algorithm do not deal with NA values: here we removed features with NA values, patients with all NA values are not removed
  view_data_new <- lapply(names(view_data), function(x) {
    tmp_data <- view_data[[x]]
    if (all(is.na(apply(tmp_data, 1, sum)))) {
      NA_sum <- apply(tmp_data, 2, sum)
      tmp_data <- tmp_data[, !is.na(NA_sum)]
    }
    return(tmp_data)
  })
  names(view_data_new) <- names(view_data)

  state <- opt_model_cancer_view_spec
  prediction_X <- view_data_new

  # standardize
  if (standardize_any == TRUE) {
    for (m in 1:P) {

      # Check features availability
      keep_pos <- stats::na.omit(match(colnames(prediction_X[[m]]), rownames(opt_xtrain_stats_cancer_view_spec[[m]]$mean)))
      keep_names <- intersect(colnames(prediction_X[[m]]), rownames(opt_xtrain_stats_cancer_view_spec[[m]]$mean))
      prediction_X[[m]] <- prediction_X[[m]][, keep_names]

      # Normalization should be done taking into account the train set
      opt_xtrain_stats_cancer_view_spec[[m]]$mean <- opt_xtrain_stats_cancer_view_spec[[m]]$mean[keep_pos, ]
      opt_xtrain_stats_cancer_view_spec[[m]]$sd <- opt_xtrain_stats_cancer_view_spec[[m]]$sd[keep_pos, ]

      prediction_X_norm <- lapply(1:K, function(k) {
        prediction_X[[m]] <- calc_z_score(
          X = prediction_X[[m]], mean = opt_xtrain_stats_cancer_view_spec[[m]]$mean[, k],
          sd = opt_xtrain_stats_cancer_view_spec[[m]]$sd[, k]
        )
      })
    }
  }

  # perform prediction
  prediction_cv <- lapply(1:K, function(k) {
    coef_matrix <- sapply(tasks, function(task) {
      state[[view_name]][[task]][, k]
    })
    rmtlr_test(prediction_X_norm[[k]], coef_matrix)
  })

  # save predictions
  prediction_cv <- lapply(1:K, function(k) {
    coef_matrix <- sapply(tasks, function(task) {
      state[[view_name]][[task]][, k]
    })
    rmtlr_test(prediction_X_norm[[k]], coef_matrix)
  })

  predictions_all_tasks_cv <- lapply(tasks, function(task) {
    prediction_task_cv <- sapply(1:K, function(k) {
      prediction_cv[[k]][, task]
    })
    return(prediction_task_cv)
  })
  names(predictions_all_tasks_cv) <- tasks
  return(predictions_all_tasks_cv)
}
