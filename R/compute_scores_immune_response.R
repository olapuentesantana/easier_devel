#' Compute published scores of immune response
#'
#' Calculates the transcriptomics-based scores
#' of hallmarks of anti-cancer immune response.
#'
#' @references Rooney, Michael S., Sachet A. Shukla, Catherine J. Wu, Gad Getz,
#' and Nir Hacohen. 2015. “Molecular and Genetic Properties of Tumors Associated
#' with Local Immune Cytolytic Activity.” Cell 160 (1): 48–61.
#' https://doi.org/10.1016/j.cell.2014.12.033.
#'
#' @references Cabrita, Rita, Martin Lauss, Adriana Sanna, Marco Donia, Mathilde
#' Skaarup Larsen, Shamik Mitra, Iva Johansson, et al. 2020. “Tertiary Lymphoid
#' Structures Improve Immunotherapy and Survival in Melanoma.”
#' Nature 577 (7791):561–65. https://doi.org/10.1038/s41586-019-1914-8.
#'
#' @references McClanahan, Mark Ayers AND Jared Lunceford AND Michael Nebozhyn
#' AND Erin Murphy AND Andrey Loboda AND David R. Kaufman AND Andrew Albright
#' AND Jonathan D. Cheng AND S. Peter Kang AND Veena Shankaran AND Sarina A.
#' Piha-Paul AND Jennifer Yearley AND Tanguy Y. Seiwert AND Antoni Ribas AND
#' Terrill K. 2017. “IFN-y–Related mRNA Profile Predicts Clinical Response
#' to PD-1 Blockade.” The Journal of Clinical Investigation 127 (8): 2930–40.
#' https://doi.org/10.1172/JCI91190.
#'
#' @references Roh, Whijae, Pei-Ling Chen, Alexandre Reuben, Christine N.
#' Spencer, Peter A. Prieto, John P. Miller, Vancheswaran Gopalakrishnan,
#' et al. 2017. “Integrated Molecular Analysis of Tumor Biopsies on Sequential
#' CTLA-4 and PD-1 Blockade Reveals Markers of Response and Resistance.”
#' Science Translational Medicine 9 (379).
#' https://doi.org/10.1126/scitranslmed.aah3560.
#'
#' @references Davoli, Teresa, Hajime Uno, Eric C. Wooten, and Stephen
#' J. Elledge. 2017. “Tumor Aneuploidy Correlates with Markers of Immune
#' Evasion and with Reduced Response to Immunotherapy.” Science 355 (6322).
#' https://doi.org/10.1126/science.aaf8399.
#'
#' @references Messina, Jane L., David A. Fenstermacher, Steven Eschrich,
#' Xiaotao Qu, Anders E. Berglund, Mark C. Lloyd, Michael J. Schell,
#' Vernon K. Sondak, Jeffrey S. Weber, and James J. Mule. 2012. “12-Chemokine
#' Gene Signature Identifies Lymph Node-Like Structures in Melanoma: Potential
#' for Patient Selection for Immunotherapy?” Scientific Reports 2 (1): 765.
#' https://doi.org/10.1038/srep00765.
#'
#' @references Auslander, Noam, Gao Zhang, Joo Sang Lee, Dennie T. Frederick,
#' Benchun Miao, Tabea Moll, Tian Tian, et al. 2018. “Robust Prediction of
#' Response to Immune Checkpoint Blockade Therapy in Metastatic Melanoma.”
#' Nature Medicine 24(10): 1545–49. https://doi.org/10.1038/s41591-018-0157-9.
#'
#' @references Fu, Yelin, Lishuang Qi, Wenbing Guo, Liangliang Jin, Kai Song,
#' Tianyi You, Shuobo Zhang, Yunyan Gu, Wenyuan Zha, and Zheng Guo. 2019. “A
#' Qualitative Transcriptional Signature for Predicting Microsatellite
#' Instability Status of Right-Sided Colon Cancer.”
#' BMC Genomics 20 (1): 769. https://doi.org/10.1186/s12864-019-6129-8.
#'
#' @references Jerby-Arnon, Livnat, Parin Shah, Michael S. Cuoco, Christopher
#' Rodman, Mei-Ju Su, Johannes C. Melms, Rachel Leeso, et al. 2018. “A Cancer
#' Cell Program Promotes t Cell Exclusion and Resistance to Checkpoint
#' Blockade.” Cell 175 (4): 984–997.e24.
#' https://doi.org/10.1016/j.cell.2018.09.006.
#'
#' @importFrom easierData get_scores_signature_genes
#'
#' @param RNA_tpm data.frame containing TPM values with
#' HGNC symbols in rows and samples in columns.
#' @param selected_scores character string with names of
#' scores of immune response to be computed. Default
#' scores are computed for: "CYT", "Roh_IS", "chemokines",
#' "Davoli_IS", "IFNy", "Ayers_expIS", "Tcell_inflamed",
#' "RIR", "TLS".
#' @param verbose logical variable indicating whether to
#' display informative messages.
#'
#' @return A numeric matrix with samples in rows and published
#' scores (gold standards) in columns.
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
#' cancer_type <- metadata(dataset_mariathasan)$cancertype
#'
#' # Select a subset of patients to reduce vignette building time.
#' pat_subset <- c(
#'   "SAM76a431ba6ce1", "SAMd3bd67996035", "SAMd3601288319e",
#'   "SAMba1a34b5a060", "SAM18a4dabbc557"
#' )
#' RNA_tpm <- RNA_tpm[, colnames(RNA_tpm) %in% pat_subset]
#'
#' # Computation of different hallmarks of anti-cancer immune responses
#' hallmarks_of_immune_response <- c(
#'   "CYT", "Roh_IS", "chemokines", "Davoli_IS", "IFNy"
#' )
#'
#' scores_immune_response <- compute_scores_immune_response(
#'   RNA_tpm = RNA_tpm,
#'   selected_scores = hallmarks_of_immune_response
#' )
compute_scores_immune_response <- function(RNA_tpm = NULL,
                                           selected_scores = c(
                                             "CYT",
                                             "Roh_IS",
                                             "chemokines",
                                             "Davoli_IS",
                                             "IFNy",
                                             "Ayers_expIS",
                                             "Tcell_inflamed",
                                             "RIR",
                                             "TLS"
                                           ),
                                           verbose = TRUE) {

  # Some checks
  if (is.null(RNA_tpm)) stop("TPM gene expression data not found")

  # Retrieve internal data
  easier_sigs <- suppressMessages(easierData::get_scores_signature_genes())

  # Check for which selected signatures appropriate functions exist
  sigs <- names(easier_sigs) %in% selected_scores
  if (verbose) {
    message(c(
      "Following scores can be computed: \n",
      paste(names(easier_sigs)[sigs], collapse = ", ")
    ), "\n")
  }

  result <- lapply(names(easier_sigs)[sigs], function(sig) {
    tryCatch(
      {
        message(sig, " score: \n")

        # Start gene re-annotation #

        ## get signature genes
        easier_spec_sig <- unique(unlist(easier_sigs[[sig]]))

        ## check for T cell inflamed
        if (sig == "Tcell_inflamed") easier_spec_sig <- easier_spec_sig[-grep("0", easier_spec_sig)]

        ## match
        literature_matches <- match(easier_spec_sig, rownames(RNA_tpm))

        if (anyNA(literature_matches)) {

          message(c(
            "Differenty named or missing signature genes : \n",
            paste(easier_spec_sig[is.na(literature_matches)], collapse = "\n")
          ))

          genes_to_annot <- easier_spec_sig[!easier_spec_sig %in% rownames(RNA_tpm)]
          out_annot <- reannotate_genes(cur_genes = genes_to_annot)

          # re-annotate gene signatures
          if (is.null(out_annot) == FALSE){

            if (length(names(easier_sigs[[sig]])) == 0){

              tmp <- easier_sigs[[sig]]
              if (any(!is.na(match(out_annot$old_names, tmp)))) {
                which_genes_old <- stats::na.omit(match(out_annot$old_names, tmp))
                which_genes_new <- !is.na(match(out_annot$old_names, tmp))
                tmp[which_genes_old] <- out_annot$new_names[which_genes_new]
                easier_sigs[[sig]] <- tmp
              }

            }else{ # sublist

              tmp <- lapply(names(easier_sigs[[sig]]), function(X) {
                tmp <- easier_sigs[[sig]][[X]]
                if (any(!is.na(match(out_annot$old_names, tmp)))) {
                  which_genes_old <- stats::na.omit(match(out_annot$old_names, tmp))
                  which_genes_new <- !is.na(match(out_annot$old_names, tmp))
                  tmp[which_genes_old] <- out_annot$new_names[which_genes_new]
                }
                return(tmp)
              })
              names(tmp) <- names(easier_sigs[[sig]])
              easier_sigs[[sig]] <- tmp

            }
          }
          message("Perfomed gene re-annotation to match signature genes  \n")
        }

        if (sig == "Tcell_inflamed") {

          # Subset RNA_tpm
          match_genes_housekeeping <- match(
            easier_sigs$Tcell_inflamed$Housekeeping.read,
            rownames(RNA_tpm)
          )
          match_genes_predictors <- match(
            easier_sigs$Tcell_inflamed$Tcell_inflamed.read,
            rownames(RNA_tpm)
          )

          # Compute score
          do.call(
            paste0("compute_", sig),
            args = list(
              housekeeping = match_genes_housekeeping,
              predictors = match_genes_predictors,
              weights = easier_sigs$Tcell_inflamed$weights,
              RNA_tpm = RNA_tpm
            ))

        } else if (sig == "IMPRES" | sig == "MSI") {

          # Subset RNA_tpm
          match_F_1 <- match(easier_sigs[[sig]]$Gene_1, rownames(RNA_tpm))
          match_F_2 <- match(easier_sigs[[sig]]$Gene_2, rownames(RNA_tpm))

          do.call("compute_IMPRES_MSI",
                  args = list(
                    sig = sig,
                    len = length(easier_sigs[sig]$Gene_1),
                    match_F_1 = match_F_1,
                    match_F_2 = match_F_2,
                    RNA_tpm = RNA_tpm
                  ))

        } else if (sig == "RIR") {

          do.call(paste0("compute_", sig), args = list(
            RNA_tpm = RNA_tpm,
            RIR_program = easier_sigs[[sig]]
          ))

        }else {

          easier_spec_sig <- unlist(easier_sigs[[sig]])
          literature_matches <- match(easier_spec_sig, rownames(RNA_tpm))
          literature_matches <- literature_matches[!is.na(literature_matches)]

          do.call(paste0("compute_", sig),
                  args = list(
                    matches = literature_matches,
                    RNA_tpm = RNA_tpm
                  ))

        }
      },
      error = function(cond) {
        cat("The following error occurred while computing sinature of ", sig, ": \n")
        df <- data.frame(rep(NA, ncol(RNA_tpm)), row.names = colnames(RNA_tpm))
        names(df)[1] <- sig
        return(df)
      }
      # ,
      # warning = function(cond) {
      # message(paste("The following warning occurred while computing sinature of", sig, ":"))
      # message(paste(cond, collapse = "/n"))
      # df <- data.frame(rep(NA, ncol(RNA_tpm)), row.names = colnames(RNA_tpm))
      # names(df)[1] <- sig
      # return(df)
      # }
    )
  })
  return(as.data.frame(result))
}
