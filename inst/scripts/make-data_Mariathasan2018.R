# Download mariathasan processed data

# Processed data used in Mariathasan, S., Turley, S., Nickles, D. et al.
# TGFβ attenuates tumour response to PD-L1 blockade by contributing to exclusion of T cells.
# Nature 554, 544–548 (2018). https://doi.org/10.1038/nature25501

# RNA-seq processed data of Mariathasan cohort is made available via IMvigor210CoreBiologies package.
# We follow instructions provided in http://research-pub.gene.com/IMvigor210CoreBiologies/)
# to derive the data.

# creating a function to clean up the data and to retrieve the data of interest.
preprocess_mariathasan <- function(cds,
                                   verbose = TRUE) {
  library(utils)
  library(Biobase)
  library(BiocGenerics)

  # Load CountDataSet object
  utils::data(cds)
  # Retrieve counts and genes information
  gene_count <- BiocGenerics::counts(cds)
  gene_info <- Biobase::fData(cds)
  clin_info <- Biobase::pData(cds)

  # Are there any duplicated symbols?
  all_duplicates <- gene_info[duplicated(gene_info$symbol), ]
  in_duplicates <- all_duplicates[!all_duplicates$symbol %in% c(NA, ""), ]
  out_duplicates <- all_duplicates[all_duplicates$symbol %in% c(NA, ""), ]

  # NCBI source (CSNK1E =1454; PTEP2-CSNK1E = 102800317)
  gene_info[which(gene_info$entrez_id == 102800317), c("symbol", "Symbol")] <- c("TPTEP2-CSNK1E", "TPTEP2-CSNK1E") # NCBI source (CSNK1E =1454; PTEP2-CSNK1E = 102800317)

  # Remove entrez id duplicates (empty symbols)
  gene_info <- gene_info[-which(gene_info$entrez_id %in% out_duplicates$entrez_id), ]
  gene_count <- gene_count[-which(rownames(gene_count) %in% out_duplicates$entrez_id), ]

  # Remove NA symbol values
  where_NA <- is.na(gene_info$symbol)
  gene_info <- gene_info[!where_NA, ]
  gene_count <- gene_count[!where_NA, ]
  rownames(gene_count) <- gene_info$symbol # Perfect match

  # obtain tpm from counts
  gene_length <- gene_info$length
  gene_count_by_len <- sweep(gene_count, 1, gene_length / 1000, FUN = "/")
  scaling_factor <- colSums(gene_count_by_len) / 1e6
  gene_tpm <- sweep(gene_count_by_len, 2, scaling_factor, FUN = "/")

  # keep only those patients with available and unambiguous response (CR and PD)
  BOR <- clin_info[, "Best Confirmed Overall Response", drop = FALSE]
  BOR_subset <- BOR[!BOR$`Best Confirmed Overall Response` == "NE", , drop = FALSE]
  BOR_subset <- BOR_subset[BOR_subset$`Best Confirmed Overall Response` %in% c("CR", "PD"), , drop = FALSE]
  gene_tpm_subset <- as.data.frame(gene_tpm[, match(rownames(BOR_subset), colnames(gene_tpm))])
  gene_count_subset <- as.data.frame(gene_count[, match(rownames(BOR_subset), colnames(gene_count))])
  clin_info_subset <- as.data.frame(clin_info[rownames(clin_info) %in% rownames(BOR_subset), , drop = FALSE])

  data <- list(counts = gene_count_subset, tpm = gene_tpm_subset, clinical = clin_info_subset)

  if (verbose) message("Mariathasan data retrieved successfully!")
  return(data)
}

# Retrieve data from package
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("biomaRt",
                       "circlize",
                       "ComplexHeatmap",
                       "corrplot",
                       "DESeq2",
                       "dplyr",
                       "DT",
                       "edgeR",
                       "ggplot2",
                       "limma",
                       "lsmeans",
                       "reshape2",
                       "spatstat",
                       "survival",
                       "plyr"))

install.packages("path/to/IMvigor210CoreBiologies_1.0.0.tar.gz", repos = NULL)
library(IMvigor210CoreBiologies)

# retrieve data
data(cds)
mariathasan_data <- preprocess_mariathasan(cds)
rm(cds)

# retrieve clinical response
clinical_data <- mariathasan_data$clinical
patient_response <- clinical_data[, "Best Confirmed Overall Response", drop=FALSE]
patient_response[,1] <- gsub("CR","R", patient_response[,1])
patient_response[,1] <- gsub("PD","NR", patient_response[,1])
# names(patient_response) <- rownames(clinical_data)

# retrieve Tumor Mutational Burden
TumorMutationalBurden <- clinical_data[, "FMOne mutation burden per MB", drop=FALSE]
# names(TumorMutationalBurden) <- rownames(clinical_data)
# # Create a external dataset object
# setClass("DatasetExample", slots=list(name="character", counts="data.frame", tpm="data.frame", response="character", TMB="numeric", cancertype="character"))
#
# dataset_mariathasan <- new("DatasetExample",
#                            name="Mariathasan",
#                            counts=mariathasan_data$counts,
#                            tpm=mariathasan_data$tpm,
#                            response=patient_response,
#                            TMB=TumorMutationalBurden,
#                            cancertype="BLCA"
# )

# ready to construct the SummarizedExperiment object
library(SummarizedExperiment)

# colData:
# matched patient's id, response to anti-PD1 and tumor mutational burden
coldata_Mariathasan2018 <- DataFrame(
  pat_id = colnames(mariathasan_data$tpm)
)
coldata_Mariathasan2018$BOR <- patient_response$`Best Confirmed Overall Response`
coldata_Mariathasan2018$TMB <- TumorMutationalBurden$`FMOne mutation burden per MB`

# cancer type from which the bulk RNA-seq was obtained
metadata_Mariathsan2018 <- DataFrame(
  cancertype = "BLCA"
)

Mariathasan2018_PDL1_treatment <- SummarizedExperiment(
  assays = SimpleList(
    counts = mariathasan_data$counts,
    tpm = mariathasan_data$tpm
  ),
  metadata = metadata_Mariathsan2018,
  colData = coldata_Mariathasan2018
)

# save
save(Mariathasan2018_PDL1_treatment, file = "inst/extdata/Mariathasan2018_PDL1_treatment.rda")


