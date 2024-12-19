#' Gene re-annotation using HGNC symbols
#'
#' Performs gene re-annotation using curated data
#' from the HGNC.
#'
#' Source code adapted from quanTIseq helper function
#' mapGenes from quantiseqr package.
#'
#' @importFrom easierData get_HGNC_annotation
#'
#' @param cur_genes character string containing gene HGNC
#' symbols to be consider for re-annotation.
#'
#' @return A data.frame with the old gene HGNC symbol and
#' the new corresponding gene HGNC symbol.
#' 
#' @export
#'
#' @examples
#' \donttest{
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
#' # Select some genes to check possible updated gene names
#' genes_to_check <- rownames(RNA_tpm)[400:450]
#' genes_info <- reannotate_genes(cur_genes = genes_to_check)
#' }
reannotate_genes <- function(cur_genes) {

  # Some checks
  if (is.null(cur_genes)) stop("Character string with gene names not found")

  # Retrieve internal data
  HGNC <- suppressMessages(easierData::get_HGNC_annotation())

  # cur_genes <- rownames(gene_expr)
  new_genes <- rep(NA, length(cur_genes))
  new_genes2 <- rep(NA, length(cur_genes))
  ind <- match(cur_genes, HGNC$ApprovedSymbol)

  # Current symbols and withdrawn ones
  genes_ind_notNA <- which(!is.na(ind))
  for (i in genes_ind_notNA) {
    genei <- cur_genes[i]
    if (HGNC$Status[ind[i]] == "Approved") {
      new_genes[i] <- cur_genes[i]
    } else if (HGNC$Status[ind[i]] == "EntryWithdrawn") {
      next
    } else {
      W_string <- "symbolwithdrawn,see"
      new_symbol <- gsub(W_string, "", HGNC$ApprovedName[ind[i]])
      new_genes2[i] <- new_symbol
    }
  }

  # Not found as symbols
  genes_ind_NA <- which(is.na(ind))
  for (i in genes_ind_NA) {
    genei <- cur_genes[i]

    # Previos symbol?
    ind1 <- grep(genei, HGNC$PreviousSymbols)
    for (i1 in ind1) {
      array1 <- unlist(strsplit(as.character(HGNC$PreviousSymbols[i1]), "|",
        fixed = TRUE
      ))
      flag1 <- length(which(array1 == genei)) > 0
      if (flag1) {
        new_symbol <- as.character(HGNC$ApprovedSymbol[i1])
        new_genes2[i] <- new_symbol
      }
    }
    # Synonym?
    ind2 <- grep(genei, HGNC$Synonyms)
    for (i2 in ind2) {
      array2 <- unlist(strsplit(as.character(HGNC$Synonyms[i2]), "|"))
      flag2 <- length(which(array2 == genei)) > 0
      if (flag2) {
        new_symbol <- as.character(HGNC$ApprovedSymbol[i2])
        new_genes2[i] <- new_symbol
      }
    }
  }

  new_genes2[which(new_genes2 %in% setdiff(new_genes, NA))] <- NA
  ind <- intersect(
    which(is.na(new_genes)),
    which(!is.na(new_genes2))
  )
  new_genes[ind] <- new_genes2[ind]

  out <- data.frame(old_names = cur_genes, new_names = new_genes)

  return(out)
}
