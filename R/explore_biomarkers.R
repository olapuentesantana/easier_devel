#' Visualization of stunning biomarkers
#'
#' This function provides an overview of relevant computed features (biomarkers), comparing responders and non-responders if known.
#' Information about the features contribution to the optimal models is also included.
#'
#' @importFrom stats aggregate
#' @importFrom reshape2 melt
#' @import rstatix
#' @importFrom ggrepel geom_text_repel
#' @import ggplot2
#' @importFrom grid unit.pmax grid.newpage grid.draw
#' @importFrom easierData get_opt_models get_intercell_networks
#'
#' @param pathways numeric matrix with pathways activity (rows = samples; columns = pathways).
#' @param immunecells numeric matrix with immune cell quantification (rows = samples; columns = cell types).
#' @param tfs numeric matrix with transcription factors activity (rows = samples; columns = transcription factors).
#' @param lrpairs numeric matrix with ligand-receptor weights (rows = samples; columns = ligand-receptor pairs).
#' @param ccpairs numeric matrix with cell-cell scores (rows = samples; columns = cell-cell pairs).
#' @param cancer_type character string indicating which cancer-specific model should be used to compute the predictions.
#' @param patient_response character vector with two factors (Non-responders = NR, Responders = R).
#' @param verbose logical flag indicating whether to display messages about the process.
#'
#' @return \itemize{
#' \item{Volcano plot displaying relevant biomarkers differentiating responders vs non-responders patients.}
#' \item{A combined plot for each type of quantitative descriptors, showing the original distribution of the features and the importance of these features for the trained models}
#' }
#'
#' @export
#'
#' @examples
#' # Load exemplary dataset (Mariathasan et al., Nature, 2018) from ExperimentHub easierData.
#' # Original processed data is available from IMvigor210CoreBiologies package.
#' library("easierData")
#' dataset_mariathasan <- easierData::get_Mariathasan2018_PDL1_treatment()
#' RNA_tpm <- dataset_mariathasan@assays@data@listData[["tpm"]]
#' RNA_counts <- dataset_mariathasan@assays@data@listData[["counts"]]
#' cancer_type <- dataset_mariathasan@metadata$cancertype
#'
#' # Computation of cell fractions
#' cell_fractions <- compute_cell_fractions(RNA_tpm = RNA_tpm)
#'
#' # Computation of pathway scores
#' pathway_activity <- compute_pathway_activity(
#'   RNA_counts = RNA_counts,
#'   remove_sig_genes_immune_response = TRUE
#' )
#'
#' # Computation of TF activity
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
#' # retrieve clinical response
#' patient_ICBresponse <- dataset_mariathasan@colData$BOR
#' names(patient_ICBresponse) <- dataset_mariathasan@colData$pat_id
#'
#' # Investigate possible biomarkers
#' explore_biomarkers(
#'   pathways = pathway_activity,
#'   immunecells = cell_fractions,
#'   lrpairs = lrpair_weights,
#'   tfs = tf_activity,
#'   ccpairs = ccpair_scores,
#'   cancer_type = cancer_type,
#'   patient_response = patient_ICBresponse
#' )
explore_biomarkers <- function(pathways = NULL,
                               immunecells = NULL,
                               tfs = NULL,
                               lrpairs = NULL,
                               ccpairs = NULL,
                               cancer_type,
                               patient_response,
                               verbose = TRUE) {
  if (missing(cancer_type)) stop("cancer type needs to be specified")
  if (missing(patient_response)) stop("patient response needs to be specified")
  if (all(is.null(pathways), is.null(immunecells), is.null(tfs), is.null(lrpairs), is.null(ccpairs))) stop("none signature specified")
  # Retrieve internal data
  opt_models <- suppressMessages(easierData::get_opt_models())
  intercell_networks <- suppressMessages(easierData::get_intercell_networks())
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

  get_biomarkers_features <- function(view, cancer_type, verbose = TRUE) {
    view_info <- view_combinations[[view]]
    view_name <- paste(names(view_info), collapse = "_")
    if (verbose) message("Examining ", view_name, " biomarkers \n")
    # ---------- #
    # Features #
    # ---------- #
    features <- as.matrix(get(view_name))
    features_z <- calc_z_score(features)
    patients <- intersect(names(patient_response), rownames(features))
    # add response labels
    response <- patient_response[patients]
    response_df <- data.frame(sample = names(response), label = response)
    features <- features[patients, ]
    features_z <- features_z[patients, ]
    features_df <- reshape2::melt(features)
    features_df_z <- reshape2::melt(features_z)
    features_df$value_z <- features_df_z$value
    names(features_df) <- c("sample", "feature", "value", "value_z")
    # Merge
    features_df <- merge(features_df, response_df)
    features_df$datatype <- view_name
    # ---------- #
    # Weights #
    # ---------- #
    opt_model_cancer_view_spec <- opt_models[[cancer_type]][[view_name]]

    my_coefs <- reshape2::melt(opt_model_cancer_view_spec)
    names(my_coefs) <- c("feature", "run", "estimate", "task")
    my_coefs$datatype <- view_name
    # remove intercept
    my_coefs <- subset(my_coefs, !feature %in% "(Intercept)")
    # calculate median across runs
    # my_coefs_task_median <- stats::aggregate(estimate ~ feature + task,
    #                                          FUN = "median", na.rm = TRUE, data = my_coefs
    # )
    # calculate median across runs and tasks
    my_coefs_median <- stats::aggregate(estimate ~ feature + datatype,
      FUN = "median", na.rm = TRUE, data = my_coefs
    )

    return(list(weights = my_coefs_median, features = features_df))
  }

  plot_list <- list()
  comparison <- do.call(rbind, lapply(1:length(view_combinations), function(ii) {
    biomarkers_weights_features <- get_biomarkers_features(ii, cancer_type)
    biomarkers_weights <- biomarkers_weights_features$weights
    biomarkers_weights$feature <- droplevels(biomarkers_weights$feature)

    # change names for cell-cell pairs
    if (unique(biomarkers_weights$datatype) == "ccpairs") {
      new_variables_cc <- c("CD8", "M", "B", "DC", "Endo", "Mast", "Fib", "Adip", "CD4", "NK", "Neu", "Mono", "Cancer")
      old_variables_cc <- c(
        "CD8+T-Cell", "Macrophages", "B-Cell", "DendriticCells", "Endothelialcells", "Mastcells", "Fibroblasts", "Adipocytes",
        "CD4+T-Cell", "NKcells", "Neutrophils", "Monocytes", "Cancercell"
      )
      tmp <- as.character(biomarkers_weights$feature)
      for (X in 1:length(new_variables_cc)) {
        tmp <- gsub(old_variables_cc[X], new_variables_cc[X], tmp, fixed = T)
      }
      biomarkers_weights$feature <- tmp
    }
    colnames(biomarkers_weights) <- c("variable", "datatype", "weight")

    biomarkers_weights_sort <- biomarkers_weights[order(abs(biomarkers_weights$weight), decreasing = TRUE), ]

    if (nrow(biomarkers_weights_sort) > 15) {
      biomarkers_weights_sort <- biomarkers_weights_sort[1:15, ]
    }
    biomarkers_weights_sort <- as.data.frame(biomarkers_weights_sort)

    features <- biomarkers_weights_features$features
    features <- features[!is.na(features$value), ]
    features$feature <- droplevels(features$feature)

    if (unique(features$datatype) == "ccpairs") {
      tmp <- features$feature
      tmp_2 <- levels(features$feature)
      for (X in 1:length(new_variables_cc)) {
        tmp <- gsub(old_variables_cc[X], new_variables_cc[X], tmp, fixed = T)
        tmp_2 <- gsub(old_variables_cc[X], new_variables_cc[X], tmp_2, fixed = T)
      }
      features$feature <- factor(tmp, levels = tmp_2)
    }

    # weights
    biomarkers_weights_sort$variable <- factor(biomarkers_weights_sort$variable, levels = unique(biomarkers_weights_sort$variable))

    biomarkers_weights_sort$cor <- sign(biomarkers_weights_sort$weight)
    biomarkers_weights_sort$cor <- gsub("-1", "-", biomarkers_weights_sort$cor, fixed = TRUE)
    biomarkers_weights_sort$cor <- gsub("1", "+", biomarkers_weights_sort$cor, fixed = TRUE)
    biomarkers_weights_sort$cor <- factor(biomarkers_weights_sort$cor, levels = unique(biomarkers_weights_sort$cor))

    # BARPLOT #
    barplot <- ggplot2::ggplot(biomarkers_weights_sort, ggplot2::aes(x = .data$variable, y = abs(.data$weight), fill = .data$cor)) +
      ggplot2::geom_bar(stat = "identity", color = "white") +
      ggplot2::scale_fill_manual(
        name = "Association sign",
        labels = c("+", "-", "0"),
        values = c("+" = "#4477AA", "-" = "#BB4444", "0" = "gray")
      ) +
      ggplot2::theme(panel.grid = ggplot2::element_blank()) +
      ggplot2::theme(
        axis.text.y = ggplot2::element_text(size = 12, color = "black"), axis.title.x = ggplot2::element_blank(), axis.title.y = ggplot2::element_text(size = 12),
        axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_line(size = 0.5, color = "black"),
        legend.position = "top", legend.direction = "horizontal",
        legend.box.background = ggplot2::element_rect(color = "black", size = 0.3),
        legend.box.margin = ggplot2::margin(0.5, 0.5, 0.5, 0.5),
        legend.text = ggplot2::element_text(size = 12),
        legend.title = ggplot2::element_text(size = 12, face = "bold", vjust = 0.5),
        panel.border = ggplot2::element_blank(), panel.background = ggplot2::element_blank(),
        plot.margin = ggplot2::unit(c(0, 0, 0.2, 0.2), "cm"), axis.line.y = ggplot2::element_line(colour = "black")
      ) +
      ggplot2::labs(y = "Biomarker weight") +
      ggplot2::labs(title = paste0(" Quantitative descriptor: ", unique(biomarkers_weights_sort$datatype)))


    features_boxplot <- subset(features, feature %in% unique(biomarkers_weights_sort$variable))
    features_boxplot$feature <- factor(as.character(features_boxplot$feature), levels = unique(biomarkers_weights_sort$variable))
    features_boxplot$label <- factor(features_boxplot$label, levels = c("NR", "R"))

    # BOXPLOT #
    boxplot <- ggplot2::ggplot(features_boxplot, ggplot2::aes(x = .data$feature, y = .data$value_z, fill = .data$label, color = .data$label)) +
      ggplot2::geom_boxplot(alpha = 0.8, outlier.shape = NA) +
      ggplot2::geom_point(position = ggplot2::position_jitterdodge(), size = 0.05) +
      ggplot2::scale_fill_manual(
        name = "Label",
        labels = levels(features_boxplot$label),
        values = c("darkgrey", "black")
      ) +
      ggplot2::scale_color_manual(
        name = "Label",
        labels = levels(features_boxplot$label),
        values = c("darkgrey", "black")
      ) +
      ggplot2::ylim(c(-5, 5)) +
      ggplot2::theme_minimal() +
      ggplot2::theme(panel.grid = ggplot2::element_blank()) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(size = 12, angle = 45, hjust = 1, color = "black"), axis.text.y = ggplot2::element_text(size = 12, color = "black"),
        axis.title.x = ggplot2::element_blank(), axis.title.y = ggplot2::element_text(size = 12), axis.ticks.x = ggplot2::element_blank(),
        legend.position = "bottom", legend.direction = "horizontal", axis.ticks.y = ggplot2::element_line(size = 0.5, color = "black"),
        legend.box.background = ggplot2::element_rect(color = "black", size = 0.3),
        legend.box.margin = ggplot2::margin(0.5, 0.5, 0.5, 0.5),
        legend.text = ggplot2::element_text(size = 12),
        legend.title = ggplot2::element_text(size = 12, face = "bold", vjust = 0.5),
        plot.margin = ggplot2::unit(c(0.2, 0, 0, 0.2), "cm"), axis.line.y = ggplot2::element_line(colour = "black")
      ) +
      ggplot2::labs(y = "Z-score")

    if (length(unique(biomarkers_weights$feature)) > 30) {
      boxplot <- boxplot + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 10))
    } else {
      boxplot <- boxplot + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 12))
    }
    # Combine plots
    g1 <- ggplot2::ggplotGrob(boxplot)
    g2 <- ggplot2::ggplotGrob(barplot)
    g <- rbind(g2, g1, size = "first")
    g$widths <- grid::unit.pmax(g1$widths, g2$widths)

    grid::grid.newpage()
    suppressWarnings(grid::grid.draw(g))
    plot_list[[ii]] <- recordPlot()

    features_names <- levels(features$feature)
    datatype_comparison <- do.call(rbind, lapply(features_names, function(x) {
      # consider only the data for that variable and separate R and NR
      tmp <- subset(features, feature == x)
      # compute average
      tmp_mean <- tapply(tmp$value, tmp$label, mean)

      if (length(unique(tmp$value)) != 1) {
        # compute wilcoxon sum rank test, effect size and sign
        stattest <- rstatix::wilcox_test(data = tmp, value ~ label)
        effsize <- rstatix::wilcox_effsize(data = tmp, value ~ label)
        sign <- sign(tmp_mean["R"] - tmp_mean["NR"])
        p_val <- stattest$p
        eff_size <- effsize$effsize
      } else {
        p_val <- 1
        eff_size <- 0
        sign <- 0
      }
      tmp_df <- data.frame(
        datatype = unique(features$datatype),
        variable = x,
        p_val = p_val,
        eff_size = eff_size,
        sign = sign
      )
      return(tmp_df)
    }))
    datatype_comparison <- merge(datatype_comparison, biomarkers_weights)
    datatype_comparison$istop <- FALSE
    datatype_comparison$istop[datatype_comparison$variable %in% biomarkers_weights_sort$variable] <- TRUE
    return(datatype_comparison)
  }))

  comparison$signedEffect <- comparison$eff_size * comparison$sign
  comparison$threshold <- as.numeric(as.factor(comparison$p_val <= 0.05))
  comparison$threshold <- factor(comparison$threshold, levels = c(1, 2), labels = c("notSign", "Sign"))

  # Add arrow in lrpairs and ccpairs
  if (any(comparison$datatype %in% c("lrpairs", "ccpairs"))) {
    tmp <- sapply(strsplit(as.character(comparison$variable)[which(comparison$datatype %in% c("lrpairs", "ccpairs"))], split = "_"), function(X) {
      return(X[1:8])
    })

    # LR pairs network
    intercell_network <- intercell_networks[["pancan"]]
    LR_pairs <- unique(paste0(intercell_network$ligands, "_", intercell_network$receptors))

    new_name <- do.call(c, lapply(1:ncol(tmp), function(X) {
      tmp_2 <- tmp[!(is.na(tmp[, X])), X]
      if (length(tmp_2) > 2) {
        pos_comb <- combn(length(tmp_2), 2)
        search <- sapply(1:ncol(pos_comb), function(X) {
          paste(tmp_2[pos_comb[, X]], collapse = "_")
        })
        keep <- search[search %in% LR_pairs]
        maj <- names(which(table(unlist(strsplit(keep, split = "_"))) > 1))
        other <- paste(tmp_2[!tmp_2 %in% maj])

        if (match(maj, tmp_2) == length(tmp_2)) {
          new_name <- paste0(paste(other, collapse = "_"), "->", maj)
        } else {
          new_name <- paste0(maj, "->", paste(other, collapse = "_"))
        }
      } else if (length(tmp_2) <= 2) {
        new_name <- paste(tmp_2, collapse = "->")
      }
      return(new_name)
    }))
    comparison$variable[which(comparison$datatype %in% c("lrpairs", "ccpairs"))] <- new_name
  }

  xminmax <- max(abs(comparison$signedEffect))
  xminmax <- xminmax + xminmax * 0.01
  ymax <- max(-log10(comparison$p_val))
  ymax <- ymax + ymax * 0.01

  volcano_plot <- ggplot2::ggplot(data = comparison, ggplot2::aes(x = .data$signedEffect, y = -log10(.data$p_val), color = .data$threshold, size = abs(.data$weight))) +
    ggplot2::geom_point(alpha = 1, ggplot2::aes(shape = as.factor(sign(weight)))) +
    ggplot2::xlim(c(-xminmax, xminmax)) +
    ggplot2::ylim(c(0, ymax)) +
    ggplot2::xlab("higher in NR          effect size          higher in R") +
    ggplot2::ylab("-log10 p-value") +
    ggplot2::ggtitle("") +
    ggplot2::scale_color_manual(values = c("notSign" = "#a6a6a6", "Sign" = "#4BA8D7"), name = "R vs NR significance") +
    ggplot2::scale_shape_manual(values = c(15, 16, 17), name = "Association sign") +
    ggplot2::scale_size_continuous(name = "Estimated weight") +
    ggplot2::theme_bw() +
    ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "longdash", colour = "#9e9e9e") +
    ggplot2::geom_vline(xintercept = 0, linetype = "solid", colour = "#9e9e9e") +
    ggplot2::theme(axis.text = ggplot2::element_text(color = "black"), axis.ticks = ggplot2::element_line(color = "black")) +
    ggplot2::theme(legend.position = "right") +
    ggrepel::geom_text_repel(
      data = subset(comparison, (threshold != "notSign")), ggplot2::aes(x = .data$signedEffect, y = -log10(.data$p_val), label = .data$variable, size = .05),
      show.legend = NA, inherit.aes = FALSE, max.overlaps = 20
    ) +
    ggplot2::labs(title = paste0(" All quantitative descriptor at once"))

  plot_list[[length(view_combinations)+1]] <- suppressWarnings(print(volcano_plot))
  return(plot_list)
}
