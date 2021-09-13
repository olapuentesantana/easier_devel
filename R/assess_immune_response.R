#' Evaluate the predictive performance of the computed immune response score
#'
#' This function generates a roc curve plot and a barplot showing the average
#' (across tasks) area under the ROC curve (AUC) values for each quantitative descriptor on
#' the patients' response provided. An ensemble model is built based on the average of the individual models.
#' Additionally, the average of the gold standard scores is used for comparison.
#'
#' @importFrom ROCR prediction performance plot
#' @importFrom grDevices recordPlot
#' @importFrom stats aggregate median sd
#' @importFrom graphics legend par title abline lines plot.new
#' @import ggplot2
#'
#' @export
#'
#' @param predictions_immune_response list containing the predictions for each quantitative descriptor and for each task.
#' @param patient_response character vector with two factors (Non-responders = NR, Responders = R).
#' @param RNA_tpm numeric matrix of patients' gene expression data as tpm values.
#' @param TMB_values numeric vector containing patients' tumor mutational burden (TMB) values.
#' @param easier_with_TMB logical flag indicating whether to apply refined approach using the combination of easier
#' predictions and tumor mutational burden.
#' @param verbose logical flag indicating whether to display messages about the process.
#'
#' @return If easier_with_TMB is set to FALSE, two figures (roc curve and bar plots) are directly saved in the path specified
#' in output_file_path. If easier_with_TMB is set to TRUE, an additional plot is returned displaying an integrated approach
#' that uses both immune response and tumor mutational burden (TMB) to predict patients' response.
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
#' # Computation of cell fractions  (Finotello et al., Genome Med, 2019)
#' cell_fractions <- compute_cell_fractions(RNA_tpm = RNA_tpm)
#'
#' # Computation of pathway scores (Holland et al., BBAGRM, 2019; Schubert et al., Nat Commun, 2018)
#' pathway_activities <- compute_pathway_activity(
#'   RNA_counts = RNA_counts,
#'   remove_sig_genes_immune_response = TRUE
#' )
#'
#' # Computation of TF activity (Garcia-Alonso et al., Genome Res, 2019)
#' tf_activities <- compute_TF_activity(
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
#'   pathways = pathway_activities,
#'   immunecells = cell_fractions,
#'   tfs = tf_activities,
#'   lrpairs = lrpair_weights,
#'   ccpairs = ccpair_scores,
#'   cancer_type = cancer_type,
#'   verbose = TRUE
#' )
#'
#' # retrieve clinical response
#' patient_ICBresponse <- dataset_mariathasan@colData$BOR
#' names(patient_ICBresponse) <- dataset_mariathasan@colData$pat_id
#'
#' # retrieve TMB
#' TMB <- dataset_mariathasan@colData$TMB
#' names(TMB) <- dataset_mariathasan@colData$pat_id
#'
#' # Assess patient-specific likelihood of response to ICB therapy
#' assess_immune_response(
#'   predictions_immune_response = predictions_immune_response,
#'   patient_response = patient_ICBresponse,
#'   RNA_tpm = RNA_tpm,
#'   TMB_values = TMB,
#'   easier_with_TMB = TRUE
#' )
assess_immune_response <- function(predictions_immune_response = NULL,
                                   patient_response = NULL,
                                   RNA_tpm,
                                   TMB_values,
                                   easier_with_TMB = FALSE,
                                   verbose = TRUE) {
  if (is.null(predictions_immune_response)) stop("None predictions found")
  if (missing(TMB_values)) {
    TMB_available <- FALSE
    easier_with_TMB <- FALSE
  } else if (easier_with_TMB == FALSE){
    TMB_available <- FALSE
  } else {
    TMB_available <- TRUE
    if (anyNA(TMB_values)) warning("NA values were found in TMB data, patients with NA values are removed from the analysis")
    message(paste0("\nConsidering ", length(TMB_values[!is.na(TMB_values)]), " patients out of ", length(TMB_values), " with available TMB"))
    patients_to_keep <- names(TMB_values[!is.na(TMB_values)])
    if (is.numeric(TMB_values)) warning("Converting TMB values into numeric")
    TMB_values <- as.numeric(TMB_values[patients_to_keep])
    names(TMB_values) <- patients_to_keep
    if (is.null(patient_response) == FALSE) patient_response <- patient_response[patients_to_keep]
    RNA_tpm <- RNA_tpm[, patients_to_keep]
  }
  # Initialize variables
  views <- names(predictions_immune_response)
  tasks <- names(predictions_immune_response[[1]])
  # All views
  all_color_views <- vector("character", length = length(views))
  all_color_views <- c(
    "#52ac68", "#6a70d7", "#bbb442",
    "#5b3788", "#72a646"
  )
  names(all_color_views) <- views
  all_color_views <- all_color_views[views]

  # ---------------------------#
  # Arrange ensemble view predictor #
  # ---------------------------#
  ensemble_df <- lapply(views, function(spec_view) {
    ensemble_df <- sapply(tasks, function(spec_task) {
      df <- predictions_immune_response[[spec_view]][[spec_task]]
      if (TMB_available) df <- df[patients_to_keep, ]
      df_runs <- rowMeans(df)
    })
    return(ensemble_df)
  })
  names(ensemble_df) <- views
  overall_df <- sapply(tasks, function(spec_task) {
    overall_df <- apply(cbind(
      ensemble_df$pathways[, spec_task],
      ensemble_df$immunecells[, spec_task],
      ensemble_df$tfs[, spec_task],
      ensemble_df$lrpairs[, spec_task],
      ensemble_df$ccpairs[, spec_task]
    ), 1, mean)
  })

  # AUC predictions (when response available)
  if (is.null(patient_response) == FALSE) {
    if (all(levels(as.factor(patient_response)) %in% c("NR", "R")) == FALSE) {
      stop("patient_response factor levels are not NR and R")
    }
    # Patient response labels
    labels <- matrix(patient_response,
      nrow = length(patient_response), ncol = length(tasks),
      dimnames = list(colnames(RNA_tpm), tasks)
    )
    # Compute scores of immune response and consider them as gold standards
    tasks_values <- compute_scores_immune_response(RNA_tpm)
    if (verbose) message("Scores of immune response computed!")
    # Assess correlation between chemokines and the other correlated tasks
    tasks_cormatrix <- cor(tasks_values)
    cor_sign <- sign(tasks_cormatrix[, "chemokines"])
    cor_sign <- cor_sign[names(cor_sign) != "chemokines"]
    if (all(cor_sign == -1)) {
      tasks_values[, "chemokines"] <- -tasks_values[, "chemokines"]
    }
    tasks_values <- as.data.frame(tasks_values)
    # ---------------------------#
    # Predictions single views #
    # ---------------------------#
    ROC_pred <- lapply(views, function(spec_view) {
      ROC_pred <- sapply(tasks, function(spec_task) {
        df <- predictions_immune_response[[spec_view]][[spec_task]]
        if (TMB_available) df <- df[patients_to_keep, ]
        # check patients match
        df <- df[match(rownames(labels), rownames(df)), ]
        df_runs <- rowMeans(df)
      })
      pred <- ROCR::prediction(ROC_pred, labels, label.ordering = c("NR", "R"))
      perf <- ROCR::performance(pred, "tpr", "fpr")
      AUC <- unlist(ROCR::performance(pred, "auc")@y.values)
      names(AUC) <- tasks
      return(list(Curve = list(perf), Barplot = list(AUC)))
    })
    names(ROC_pred) <- views
    # Collect derived signatures predictions into data.frame #
    AUC_data <- do.call(rbind, lapply(views, function(spec_view) {
      return(data.frame(
        AUC = as.numeric(unlist(ROC_pred[[spec_view]]$Barplot)),
        View = spec_view,
        Task = names(unlist(ROC_pred[[spec_view]]$Barplot)),
        Run = "average"
      ))
    }))
    # Average across tasks
    AUC_mean_sd_run_tasks <- do.call(
      data.frame,
      aggregate(AUC ~ View, data = AUC_data, FUN = function(x) c(mean = mean(x), sd = sd(x)))
    )
    # ---------------------------#
    # Predictions ensemble view #
    # ---------------------------#
    # check patients match
    overall_df <- overall_df[match(rownames(labels), rownames(overall_df)), ]
    pred <- ROCR::prediction(overall_df, labels, label.ordering = c("NR", "R"))
    perf <- ROCR::performance(pred, "tpr", "fpr")
    AUC <- unlist(ROCR::performance(pred, "auc")@y.values)
    names(AUC) <- tasks
    ensemble_ROC_pred <- list(ensemble = list(Curve = list(perf), Barplot = list(AUC)))
    # Collect predictions into data.frame #
    AUC_data_ensemble <- data.frame(
      AUC = as.numeric(unlist(ensemble_ROC_pred$ensemble$Barplot)),
      View = "ensemble",
      Task = names(unlist(ensemble_ROC_pred$ensemble$Barplot))
    )
    # Average across tasks
    AUC_mean_sd_ensemble_run_tasks <- do.call(
      data.frame,
      aggregate(AUC ~ View, data = AUC_data_ensemble, FUN = function(x) c(mean = mean(x), sd = sd(x)))
    )
    # ---------------------------#
    # Predictions gold standards #
    # ---------------------------#
    # check patients match
    tasks_values <- tasks_values[match(rownames(labels), rownames(tasks_values)), ]
    pred <- ROCR::prediction(tasks_values, labels[, colnames(tasks_values)], label.ordering = c("NR", "R"))
    perf <- ROCR::performance(pred, "tpr", "fpr")
    AUC <- unlist(ROCR::performance(pred, "auc")@y.values)
    names(AUC) <- colnames(tasks_values)
    goldstandard_ROC_pred <- list(gold_standard = list(Curve = list(perf), Barplot = list(AUC)))

    AUC_data_goldstandard <- data.frame(
      AUC = as.numeric(unlist(goldstandard_ROC_pred$gold_standard$Barplot)),
      View = "gold_standard",
      Task = names(unlist(goldstandard_ROC_pred$gold_standard$Barplot))
    )
    # Average across tasks
    AUC_mean_sd_goldstandard_run_tasks <- do.call(
      data.frame,
      aggregate(AUC ~ View, data = AUC_data_goldstandard, FUN = function(x) c(mean = mean(x), sd = sd(x)))
    )
    # ------------------------------------------------------#
    # Combine single views, ensemble and gold standard
    # ------------------------------------------------------#
    AUC_mean_sd_all_run_tasks <- rbind(AUC_mean_sd_run_tasks, AUC_mean_sd_ensemble_run_tasks, AUC_mean_sd_goldstandard_run_tasks)
    ROC_all_run_tasks <- c(ROC_pred, ensemble_ROC_pred, goldstandard_ROC_pred)
    names(ROC_all_run_tasks) <- c(names(all_color_views), "ensemble", "gold_standard")
    AUC_mean_sd_all_run_tasks$View <- factor(AUC_mean_sd_all_run_tasks$View,
      levels = c(names(all_color_views), "ensemble", "gold_standard")
    )

    # Colors gold standard and ensemble
    color_gold_standard <- "gray82"
    names(color_gold_standard) <- "gold_standard"
    color_ensemble <- "gold2"
    names(color_ensemble) <- "ensemble"
    # ------------------------------------------------------#
    # Add TMB if available
    # ------------------------------------------------------#
    if (TMB_available == TRUE) {
      # check patients match
      TMB_values <- TMB_values[match(rownames(labels), names(TMB_values))]
      pred <- ROCR::prediction(TMB_values, labels[, 1], label.ordering = c("NR", "R"))
      perf <- ROCR::performance(pred, "tpr", "fpr")
      AUC <- unlist(ROCR::performance(pred, "auc")@y.values)
      TMB_ROC_pred <- list(TMB = list(Curve = list(perf), Barplot = list(AUC)))

      # Collect predictions into data.frame #
      AUC_data_TMB <- data.frame(
        AUC = as.numeric(unlist(TMB_ROC_pred$TMB$Barplot)),
        View = "TMB",
        Task = "TMB"
      )
      # Average across tasks
      AUC_mean_sd_TMB_run_tasks <- do.call(
        data.frame,
        aggregate(AUC ~ View, data = AUC_data_TMB, FUN = function(x) c(mean = mean(x), sd = sd(x)))
      )
      AUC_mean_sd_all_run_tasks <- rbind(AUC_mean_sd_all_run_tasks, AUC_mean_sd_TMB_run_tasks)
      ROC_all_run_tasks <- c(ROC_all_run_tasks, TMB_ROC_pred)
      names(ROC_all_run_tasks) <- c(names(all_color_views), "ensemble", "gold_standard", "TMB")
      AUC_mean_sd_all_run_tasks$View <- factor(AUC_mean_sd_all_run_tasks$View,
        levels = c(names(all_color_views), "ensemble", "gold_standard", "TMB")
      )
      # Colors TMB
      color_TMB <- "salmon"
      names(color_TMB) <- "salmon"
    }
    plot_list <- list()
    # *******************************************
    # Plot AUC values using barplot
    # *******************************************
    n_R <- table(patient_response)[["R"]]
    n_NR <- table(patient_response)[["NR"]]
    gg <- ggplot2::ggplot(AUC_mean_sd_all_run_tasks, ggplot2::aes(x = .data$View, y = round(.data$AUC.mean, 2), fill = .data$View)) +
      ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge(), color = "white") +
      if (TMB_available) {
        ggplot2::scale_fill_manual(values = c(
          as.vector(all_color_views), as.vector(color_ensemble),
          as.vector(color_gold_standard), as.vector(color_TMB)
        ), guide = "none")
      } else {
        ggplot2::scale_fill_manual(values = c(
          as.vector(all_color_views), as.vector(color_ensemble),
          as.vector(color_gold_standard)
        ), guide = "none")
      }
    ggg <- gg + ggplot2::scale_x_discrete(labels = c(
      "ensemble" = "Ensemble",
      "immunecells" = "Cell fractions",
      "pathways" = "Pathways",
      "tfs" = "TFs",
      "lrpairs" = "LR pairs",
      "ccpairs" = "CC pairs",
      "gold_standard" = "Tasks (gold standard)"
    )) +
      ggplot2::theme(panel.grid = ggplot2::element_blank(), panel.background = ggplot2::element_rect(fill = NA)) +
      ggplot2::theme_bw() +
      ggplot2::ylim(0, 1) +
      ggplot2::ylab("Area under the curve (AUC)") +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = round(.data$AUC.mean, 2) - .data$AUC.sd, ymax = round(.data$AUC.mean, 2) + .data$AUC.sd), width = .3, color = "black", position = ggplot2::position_dodge(0.9)) +
      ggplot2::geom_text(ggplot2::aes(label = round(.data$AUC.mean, 2)), stat = "identity", color = "black", size = 4, angle = 90, hjust = -0.5, position = ggplot2::position_dodge(0.9)) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(size = 12, angle = 45, vjust = 1, hjust = 1, color = "black"),
        axis.text.y = ggplot2::element_text(size = 12, color = "black"),
        axis.title.y = ggplot2::element_text(size = 12), axis.title.x = ggplot2::element_blank(),
        legend.position = "right", legend.text = ggplot2::element_text(size = 10),
        legend.box.background = ggplot2::element_rect(color = "black", size = 0.3),
        legend.box.margin = ggplot2::margin(0.5, 0.5, 0.5, 0.5)
      ) +
      ggplot2::labs(title = paste0(" n=", length(patient_response), " (R=", n_R, "; ", "NR=", n_NR, ")"))

    plot_list[[1]] <- print(ggg)

    # *******************************************
    # Plot ROC curve
    # *******************************************
    all_colors <- c(all_color_views, color_ensemble, color_gold_standard)
    graphics::par(cex.axis = 1.3, mar = c(5, 4, 2, 2), col.lab = "black", pty="s")

    # Single views & Gold Standard
    ROCR::plot(ROC_all_run_tasks[[1]]$Curve[[1]],
      avg = "threshold", col = all_colors[1], lwd = 2, type = "S",
      cex.lab = 1.3, ylab = "True Positive Rate", xlab = "False Positive Rate"
    )

    sapply(setdiff(names(ROC_all_run_tasks)[2:length(names(ROC_all_run_tasks))], "TMB"), function(descriptor) {
      # Single views
      ROCR::plot(ROC_all_run_tasks[[descriptor]]$Curve[[1]],
        avg = "threshold", col = all_colors[descriptor], lwd = 2, type = "S",
        cex.lab = 1.3, ylab = "True Positive Rate", xlab = "False Positive Rate", add = TRUE
      )
    })

    legend_text <- do.call(c, lapply(setdiff(names(ROC_all_run_tasks)[1:length(names(ROC_all_run_tasks))], "TMB"), function(descriptor) {
      paste0(descriptor," (AUC=", round(subset(AUC_mean_sd_all_run_tasks, View == descriptor)$AUC.mean, 2),")")
    }))

    # Ensemble along different tasks
    ROCR::plot(ROC_all_run_tasks[["ensemble"]]$Curve[[1]],
               avg = "threshold", col = all_colors["ensemble"], lwd = 2, type = "S", lty=1,
               cex.lab = 1.3, ylab = "True Positive Rate", xlab = "False Positive Rate", add=TRUE
    )

    # Gold standard along different tasks
    ROCR::plot(ROC_all_run_tasks[["gold_standard"]]$Curve[[1]],
              col = all_colors["gold_standard"], lwd = 2, type = "S", lty = 3,
               cex.lab = 1.3, ylab = "True Positive Rate", xlab = "False Positive Rate", add=TRUE
    )

    # TMB
    if (TMB_available) {
      ROCR::plot(ROC_all_run_tasks$TMB$Curve[[1]],
        col = color_TMB, lwd = 2, type = "S", cex.lab = 1.3,
        lty = 1, add = TRUE
      )

      legend_text <- do.call(c, lapply(names(ROC_all_run_tasks)[1:length(names(ROC_all_run_tasks))], function(descriptor) {
        paste(descriptor, "(", round(subset(AUC_mean_sd_all_run_tasks, View == descriptor)$AUC.mean, 2), ")")
      }))
      graphics::legend(
        x = "topleft",
        legend = legend_text,
        col = c(all_colors, color_TMB), lty = 1, lwd = 2, cex = 0.6, bty = "n"
      )
    } else {
      graphics::legend(
        x = "topleft",
        legend = legend_text,
        col = all_colors, lty = 1, lwd = 2, cex = 0.6, bty = "n"
      )
    }
    plot_list[[2]] <- grDevices::recordPlot()

    # *******************************************
    # Integrated score
    # *******************************************
    if (easier_with_TMB == TRUE) {
      rp_df <- data.frame(
        response = patient_response,
        prediction_easier = apply(overall_df, 1, mean),
        TMB = TMB_values
      )
      # Categorize TMB
      if (length(unique(rp_df$TMB)) > 3) {
        rp_df$TMBcat <- categorize_TMB(rp_df$TMB)
        # if I want specify the thresholds
        # rp.df$TMB <- categorize.TMB(rp.df$TMB, thresholds = c(100,400))
      }
      # compute the integrated score for different penalties #
      pred_combined <- rp_df$prediction_easier
      AUC_combined_v <- sapply(seq(from = 0, to = 1, by = 0.1), function(p) {
        pred_combined[rp_df$TMBcat == 1] <- pred_combined[rp_df$TMBcat == 1] - p
        pred_combined[rp_df$TMBcat == 3] <- pred_combined[rp_df$TMBcat == 3] + p

        pred <- ROCR::prediction(pred_combined, rp_df$response)
        AUC_combined <- ROCR::performance(pred, measure = "auc")
        AUC_combined_v <- AUC_combined@y.values[[1]]
      })
      # compute the integrated score as weighted average #
      linear_func <- function(x) {
        min_x <- min(x)
        max_x <- max(x)
        x01 <- (x - min_x) / (max_x - min_x)
        return(x01)
      }
      pred_lin <- linear_func(rp_df$prediction_easier)
      TMB_lin <- linear_func(rp_df$TMB)
      AUC_averaged_v <- sapply(seq(from = 0, to = 1, by = 0.1), function(p) {
        pred_averaged <- apply(cbind((1 - p) * pred_lin, (p) * TMB_lin), 1, mean)
        pred <- ROCR::prediction(pred_averaged, rp_df$response)
        AUC_averaged <- ROCR::performance(pred, measure = "auc")
        AUC_averaged_v <- AUC_averaged@y.values[[1]]
      })

      pred <- ROCR::prediction(rp_df$prediction_easier, rp_df$response)
      AUC_easier <- ROCR::performance(pred, measure = "auc")
      AUC_easier_v <- AUC_easier@y.values[[1]]

      graphics::par(cex.axis = 1.3, mar = c(5, 4, 2, 2), col.lab = "black", pty="s")
      plot(seq(from = 0, to = 1, by = 0.1), AUC_combined_v,
           xlab = "Penalty or Relative weight", ylab = "Area under the curve (AUC)",
           type = "b", col = "#c15050", lty = 1, pch = 19, lwd = 2, ylim = c(0, 1),
           cex.lab = 1.3)
      graphics::lines(seq(from = 0, to = 1, by = 0.1), AUC_averaged_v,
            xlab = "Penalty or Relative weight", ylab = "Area under the curve (AUC)",
            type = "b", col = "#693c72", lty = 1, pch = 19, lwd = 2, ylim = c(0, 1),
            cex.lab = 1.3)
      # TMB
      graphics::abline(h = AUC_mean_sd_TMB_run_tasks$AUC.mean, col = color_TMB)
      # easier (ensemble)
      graphics::abline(h = AUC_easier_v, col = color_ensemble)
      graphics::legend(
        x = "topright",
        legend = c("Penalized score", "Weighted average", "EaSIeR", "TMB"),
        col = c(
          "#c15050", "#693c72", as.vector(color_ensemble), color_TMB
        ), lty = 1, lwd = 2, cex = 0.6, bty = "n"
      )
      plot_list[[3]] <- grDevices::recordPlot()
    }
    return(plot_list)
  } else {
    # *******************************************
    # what if patients' response is not provided
    # *******************************************
    if (verbose) message("Scoring patients' as real response is not provided\n")
    plot_list <- list()

    rp_df <- data.frame(
      #response = patient_response,
      prediction_easier = apply(overall_df, 1, mean),
      patient = rownames(overall_df)
    )

    # color_response <- rp_df$response
    # names(color_response) <- rp_df$patient
    # color_response <- gsub("NR", "red", color_response)
    # color_response <- gsub("R", "blue", color_response)

    all_scores_df <- reshape2::melt(overall_df)
    # names(all_scores_df) <- c("response", "patient", "approach", "pred")
    names(all_scores_df) <- c("patient", "approach", "pred")
    # tmp_easier <- data.frame(patient = rownames(overall_df),
    #                          approach = "easier",
    #                          pred = apply(overall_df, 1, mean))
    # all_scores_df <- rbind(all_scores_df, tmp_easier)


    all_scores_df$approach <- factor(all_scores_df$approach, levels = tasks)
    #colors_approach <- c("#ff7a4a")
    # sort patients
    easier <- aggregate(pred ~ patient, data = all_scores_df, FUN = "median")
    easier <- easier[match(sort(easier$pred, decreasing = TRUE), easier$pred), ]
    all_scores_df$patient <- factor(all_scores_df$patient, levels = as.character(easier$patient))
    ordering <- all_scores_df$patient

    score_plot <- ggplot2::ggplot(all_scores_df, ggplot2::aes(x = pred, y = patient)) +
      ggplot2::geom_boxplot() +
      ggplot2::theme(panel.grid = ggplot2::element_blank(), panel.background = ggplot2::element_rect(fill = NA)) +
      ggplot2::geom_vline(xintercept = 0) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(size = 12, face = "bold", angle = 0, vjust = 0.5, hjust = 0.5, color = "black"),
        axis.title =  ggplot2::element_text(size = 14),
        axis.text.y = ggplot2::element_text(size = 12), axis.ticks.y = ggplot2::element_blank(), panel.background = ggplot2::element_rect(fill = "white"),
        legend.position = "right", legend.direction = "vertical", panel.grid.major.y = ggplot2::element_line(linetype = 1, colour = "gray"),
        legend.text = ggplot2::element_text(size = 10), legend.title = ggplot2::element_text(size = 10, face = "bold", vjust = 0.5),
        strip.background = ggplot2::element_rect(fill = "white", colour = "white"), strip.text = ggplot2::element_text(size = 10)
      ) +
      ggplot2::labs(x = "easier prediction score", y = "patients") +
      ggplot2::guides(fill = "none", color = "none", shape = "none")

    plot_list[[1]] <- print(score_plot)

    if (easier_with_TMB == TRUE) {

      # Add TMB #
      rp_df$TMB <- TMB_values

      # Categorize TMB #
      if (length(unique(rp_df$TMB)) > 3) {
        rp_df$TMBcat <- categorize_TMB(rp_df$TMB)
        # rp.df$TMB <- categorize_TMB(rp.df$TMB, thresholds = c(100,400)) # if I want specify the thresholds
      }

      # Compute the integrated score for different penalties #
      pred_combined <- rp_df$prediction_easier
      names(pred_combined) <- rp_df$patient
      pred_combined_rf <- sapply(seq(from = 0, to = 1, by = 0.1), function(p) {
        pred_combined[rp_df$TMBcat == 1] <- pred_combined[rp_df$TMBcat == 1] - p
        pred_combined[rp_df$TMBcat == 3] <- pred_combined[rp_df$TMBcat == 3] + p
        return(pred_combined)
      })

      # Compute the integrated score as weighted average #
      linear_func <- function(x) {
        min_x <- min(x)
        max_x <- max(x)
        x01 <- (x - min_x) / (max_x - min_x)
        return(x01)
      }
      pred_lin <- linear_func(rp_df$prediction_easier)
      names(pred_lin) <- rp_df$patient
      TMB_lin <- linear_func(rp_df$TMB)
      names(TMB_lin) <- rp_df$patient

      pred_averaged_rf <- sapply(seq(from = 0, to = 1, by = 0.1), function(p) {
        pred_averaged <- apply(cbind((1 - p) * pred_lin, (p) * TMB_lin), 1, mean)
      })

      rp_df$weighted_average <- pred_averaged_rf[, 6] # penalty 0.5
      rp_df$penalized_score <- NA
      rp_df$penalized_score[match(rownames(pred_combined_rf), rp_df$patient)] <- pred_combined_rf[, 6] # penalty 0.5
      rp_df$patient <- paste0(rownames(rp_df), " (TMBcat=", rp_df$TMBcat, ")")
      rp_df$TMB <- NULL
      rp_df$TMBcat <- NULL
      #rp_df$penalized_score <- NULL
      #rp_df$prediction_easier <- NULL

      all_scores_df <- reshape2::melt(rp_df)
      # names(all_scores_df) <- c("response", "patient", "approach", "pred")
      names(all_scores_df) <- c("patient", "approach", "pred")
      colors_approach <- c("#cc8100", "#00bc93", "#ff5b61")

      # sort patients
      all_scores_df$patient <- factor(all_scores_df$patient, levels = all_scores_df$patient[match(levels(ordering), sapply(strsplit(all_scores_df$patient, split = " "), head, 1))])

      # color_response <- color_response[match(sapply(strsplit(levels(all_scores_df$patient), split = " "), head, 1), names(color_response))]
      # names(color_response) <- levels(all_scores_df$patient)

      rf_score_plot <- ggplot2::ggplot(all_scores_df, ggplot2::aes(x = pred, y = patient, shape = approach)) +
        ggplot2::geom_point(ggplot2::aes(fill = approach, colour = approach), size = 3) +
        ggplot2::scale_fill_manual(
          name = "Approach",
          labels = as.character(unique(all_scores_df$approach)),
          values = colors_approach
        ) +
        ggplot2::scale_shape_manual(
          name = "Approach",
          labels = as.character(unique(all_scores_df$approach)),
          values = c(21,22,23)
        ) +
        ggplot2::scale_color_manual(
          name = "Approach",
          labels = as.character(unique(all_scores_df$approach)),
          values = colors_approach
        ) +
        ggplot2::theme(panel.grid = ggplot2::element_blank(), panel.background = ggplot2::element_rect(fill = NA)) +
        ggplot2::geom_vline(xintercept = 0) +
        #ggplot2::theme_linedraw() +
        # ggplot2::xlim(-3, 3) +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(size = 12, face = "bold", angle = 0, vjust = 0.5, hjust = 0.5, color = "black"),
          axis.title =  ggplot2::element_text(size = 14),
          axis.text.y = ggplot2::element_text(size = 12), axis.ticks.y = ggplot2::element_blank(), panel.background = ggplot2::element_rect(fill = "white"),
          legend.position = "top", legend.direction = "horizontal", panel.grid.major.y = ggplot2::element_line(linetype = 1, colour = "gray"),
          legend.text = ggplot2::element_text(size = 10), legend.title = ggplot2::element_text(size = 10, face = "bold", vjust = 0.5),
          strip.background = ggplot2::element_rect(fill = "white", colour = "white"), strip.text = ggplot2::element_text(size = 10),
          legend.justification = c("right", "top"), legend.key = element_rect(fill = "white")
        ) +
        ggplot2::labs(x = "prediction", y = "patients")

      plot_list[[2]] <- print(rf_score_plot)
    }
    return(plot_list)
  }
}
