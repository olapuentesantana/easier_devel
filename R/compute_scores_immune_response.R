#' Compute published transcriptomics-based scores of hallmarks of anti-cancer immune response
#'
#' This function computes the scores of immune response as indicated by the user.
#'
#' @importFrom easierData get_scores_signature_genes
#'
#' @param RNA_tpm data.frame containing TPM values with HGNC symbols in rows and samples in columns.
#' @param selected_scores character string with names of scores of immune response to be computed.
#' Default scores are computed for: "CYT", "Roh_IS", "chemokines", "Davoli_IS", "IFNy", "Ayers_expIS", "Tcell_inflamed", "RIR", "TLS".
#' @param verbose logical variable indicating whether to display informative messages.
#'
#' @return A numeric matrix with samples in rows and gold standard scores in columns.
#' @export
#'
#' @examples
#' # Load exemplary dataset (Mariathasan et al., Nature, 2018) from easierData.
#' # Original processed data is available from IMvigor210CoreBiologies package.
#' library("easierData")
#' dataset_mariathasan <- easierData::get_Mariathasan2018_PDL1_treatment()
#' RNA_tpm <- dataset_mariathasan@assays@data@listData[["tpm"]]
#'
#' # Computation of different hallmarks of anti-cancer immune responses
#' hallmarks_of_immune_response <- c(
#'   "CYT", "Roh_IS", "chemokines", "Davoli_IS", "IFNy",
#'   "Ayers_expIS", "Tcell_inflamed", "RIR", "TLS"
#' )
#' scores_immune_response <- compute_scores_immune_response(
#'   RNA_tpm = RNA_tpm,
#'   selected_scores = hallmarks_of_immune_response
#' )
compute_scores_immune_response <- function(RNA_tpm, selected_scores = c("CYT", "Roh_IS", "chemokines", "Davoli_IS", "IFNy", "Ayers_expIS", "Tcell_inflamed", "RIR", "TLS"), verbose = TRUE) {

  # Some checks
  if (is.null(RNA_tpm)) stop("TPM gene expression data not found")

  # Retrieve internal data
  easier_sigs <- suppressMessages(easierData::get_scores_signature_genes())

  # Check for which selected signatures appropriate functions exist
  sigs <- names(easier_sigs) %in% selected_scores
  if (verbose) message(c("Following scores can be computed: \n", paste(names(easier_sigs)[sigs], collapse = ", ")), "\n")

  result <- lapply(names(easier_sigs)[sigs], function(sig) {
    tryCatch(
      {
        if (sig == "Tcell_inflamed") {
          if (any(rownames(RNA_tpm) %in% "C14orf102")) {
            message("Gene name changed: NRDE2 is approved symbol, not C14orf102", "\n")
            rownames(RNA_tpm)[rownames(RNA_tpm) %in% "C14orf102"] <- "NRDE2"
          }

          # Subset RNA_tpm
          match_genes_housekeeping <- match(easier_sigs$Tcell_inflamed$Housekeeping.read, rownames(RNA_tpm))
          match_genes_predictors <- match(easier_sigs$Tcell_inflamed$Tcell_inflamed.read, rownames(RNA_tpm))

          if (anyNA(c(match_genes_housekeeping, match_genes_predictors))) {
            tmp <- c(easier_sigs$Tcell_inflamed$Housekeeping.read, easier_sigs$Tcell_inflamed$Tcell_inflamed.read)
            message(c(paste0("Differenty named or missing signature genes for ", sig, ": \n"), paste(tmp[!tmp %in% rownames(RNA_tpm)], collapse = "\n")), "\n")
            match_genes_housekeeping <- match_genes_housekeeping[!is.na(match_genes_housekeeping)]
            match_genes_predictors <- match_genes_predictors[!is.na(match_genes_predictors)]
          }
          do.call(
            paste0("compute_", sig),
            args = list(
              housekeeping = match_genes_housekeeping,
              predictors = match_genes_predictors,
              weights = easier_sigs$Tcell_inflamed$weights,
              RNA_tpm = RNA_tpm
            )
          )
        } else if (sig == "IMPRES" | sig == "MSI") {
          read <- unique(unlist(easier_sigs[[sig]])) # 15 genes

          # EQUIVALENT : "VISTA" = "C10orf54", "PDL-1" = "CD274", "TIM-3" = "HAVCR2",
          # "PD-1" = "PDCD1", "HVEM" = "TNFRSF14", "OX40L" = "TNFSF4", "CD137L" = "TNFSF9"

          # Some genes might have other name: case for "C10orf54", it's called "VSIR", be careful
          if (any(rownames(RNA_tpm) %in% "VSIR") & sig == "IMPRES") {
            message("Gene name changed: C10orf54 instead of VSIR", "\n")
            rownames(RNA_tpm)[rownames(RNA_tpm) %in% "VSIR"] <- "C10orf54"
          }

          # Some genes might have other name: case for "CCRN4L", it's called "NOCT", be careful
          if (any(rownames(RNA_tpm) %in% "CCRN4L") & sig == "MSI") {
            message("Gene name changed: NOCT is approved symbol, not CCRN4L", "\n")
            rownames(RNA_tpm)[rownames(RNA_tpm) %in% "CCRN4L"] <- "NOCT"
          }

          # Subset RNA_tpm
          match_F_1 <- match(easier_sigs[[sig]]$Gene_1, rownames(RNA_tpm))
          match_F_2 <- match(easier_sigs[[sig]]$Gene_2, rownames(RNA_tpm))

          if (anyNA(c(match_F_1, match_F_2))) {
            message(c("Differenty named or missing signature genes : \n", paste(read[!read %in% rownames(RNA_tpm)], collapse = "\n")))
          }

          do.call(
            "compute_IMPRES_MSI",
            args = list(
              sig = sig,
              len = length(easier_sigs[sig]$Gene_1),
              match_F_1 = match_F_1,
              match_F_2 = match_F_2,
              RNA_tpm = RNA_tpm,
              verbose = verbose
            )
          )
        } else {
          easier_spec_sig <- unlist(easier_sigs[[sig]])

          # Literature genes
          literature_matches <- match(easier_spec_sig, rownames(RNA_tpm))

          if (anyNA(literature_matches)) {
            message(c(
              paste0("Differenty named or missing signature genes for ", sig, ": \n"),
              paste(easier_spec_sig[!easier_spec_sig %in% rownames(RNA_tpm)], collapse = ", ")
            ), "\n")

            # Re-annotate genes not found
            out_annot <- reannotate_genes(easier_spec_sig[!easier_spec_sig %in% rownames(RNA_tpm)])
            easier_spec_sig[match(out_annot$old_names[!is.na(out_annot$new_names)], easier_spec_sig)] <- out_annot$new_names[!is.na(out_annot$new_names)]
            message(
              "After gene re-annotation, differently named or missing signature genes for ", sig, ": \n",
              paste(easier_spec_sig[!easier_spec_sig %in% rownames(RNA_tpm)], collapse = ", "), "\n"
            )
            literature_matches <- match(easier_spec_sig, rownames(RNA_tpm))
            literature_matches <- literature_matches[!is.na(literature_matches)]
          }

          if (sig == "RIR") {

            # Modify new RIR_sig
            tmp_RIR_sig <- sapply(names(easier_sigs[[sig]]), function(X) {
              if (any(is.na(match(out_annot$old_names, easier_sigs[[sig]][[X]])) == FALSE)) {
                easier_sigs[[sig]][[X]][stats::na.omit(match(out_annot$old_names, easier_sigs[[sig]][[X]]))] <- out_annot$new_names[!is.na(match(out_annot$old_names, easier_sigs[[sig]][[X]]))]
              }
              return(easier_sigs[[sig]][[X]])
            })
            easier_sigs[[sig]] <- tmp_RIR_sig

            do.call(paste0("compute_", sig), args = list(RNA_tpm = RNA_tpm, RIR_program = easier_sigs[[sig]]))
          } else {
            do.call(paste0("compute_", sig), args = list(matches = literature_matches, RNA_tpm = RNA_tpm))
          }
        }
      },
      error = function(cond) {
        message(paste("The following error occurred while computing sinature of", sig, ":"))
        message(paste(cond, collapse = "/n"))
        df <- data.frame(rep(NA, ncol(RNA_tpm)), row.names = colnames(RNA_tpm))
        names(df)[1] <- sig
        return(df)
      }
      #,
      #warning = function(cond) {
        #message(paste("The following warning occurred while computing sinature of", sig, ":"))
        #message(paste(cond, collapse = "/n"))
        #df <- data.frame(rep(NA, ncol(RNA_tpm)), row.names = colnames(RNA_tpm))
        #names(df)[1] <- sig
        #return(df)
      #}
    )
   }
  )
  return(as.data.frame(result))
}
