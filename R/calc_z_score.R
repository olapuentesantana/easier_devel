#' Perform Z-score normalization
#'
#' This function performs z-score normalization on a numeric matrix per column.
#' Mean and standard deviation values to be used to calculate z-score values, otherwise they are
#' calculated based on the input matrix.
#'
#' @importFrom matrixStats colSds
#'
#' @export
#'
#' @param X numeric matrix.
#' @param mean numeric vector with mean values.
#' @param sd numeric vector with standard deviation values.
#'
#' @return A numeric matrix with values as z-scores.
#'
#' @examples
#' # Load exemplary dataset (Mariathasan et al., Nature, 2018) from easierData.
#' # Original processed data is available from IMvigor210CoreBiologies package.
#' library("easierData")
#' dataset_mariathasan <- easierData::get_Mariathasan2018_PDL1_treatment()
#' RNA_tpm <- dataset_mariathasan@assays@data@listData[["tpm"]]
#'
#' # apply z-score normalization
#' tpm_zscore <- calc_z_score(t(RNA_tpm))
calc_z_score <- function(X,
                         mean,
                         sd) {
  X_scale <- matrix(0, nrow(X), ncol(X), dimnames = list(rownames(X), colnames(X)))

  if (missing(mean) & missing(sd)) {
    mean_X <- colMeans(X, na.rm = TRUE)
    sd_X <- matrixStats::colSds(as.matrix(X), na.rm = TRUE)
    X_scale <- sweep(X, 2, mean_X, FUN = "-")
    X_scale <- sweep(X_scale, 2, sd_X, FUN = "/")
  } else {
    mean <- mean[na.omit(match(colnames(X), names(mean)))]
    sd <- sd[na.omit(match(colnames(X), names(sd)))]

    X <- X[, na.omit(match(names(sd), colnames(X)))]
    X_scale <- sweep(X, 2, mean, FUN = "-")
    X_scale <- sweep(X_scale, 2, sd, FUN = "/")
  }
  return(as.matrix(X_scale))
}
