## validate with `ExperimentHubData::makeExperimentHubMetadata()`
## (do this above pkg directory)

Mariathasan2018_PDL1_treatment <- data.frame(
  stringsAsFactors = FALSE,
  Title = "Mariathasan2018_PDL1_treatment",
  Description = paste(
    "Processed RNA-seq data, quantified as counts and TPM,",
    "for bladder cancer samples pre anti-PDL1 immunotherapy.",
    "Information on patient's tumor mutational burden",
    "and response to anti-PD1 treatment is also provided.",
    "Represented as a SummarizedExperiment;",
    "derived from IMvigor210CoreBiologies R package",
    "(Mariathasan et al., Nature, 2018) under the CC-BY license."
  ),
  BiocVersion = "3.13",
  Genome = NA,
  SourceType = "tar.gz",
  SourceUrl = "http://research-pub.gene.com/IMvigor210CoreBiologies/",
  SourceVersion = "Feb 14 2018",
  Species = "Homo sapiens",
  TaxonomyId = "9606",
  Coordinate_1_based = NA,
  DataProvider = paste(
    "Mariathasan S, Turley S, Nickles D et al., “TGF-b",
    "attenuates tumor response to PD-L1 blockade by contributing to exclusion of T cells."
  ),
  Maintainer = "Federico Marini <marinif@uni-mainz.de>",
  RDataClass = "SummarizedExperiment",
  DispatchClass = "Rda",
  RDataPath = "easier/Mariathasan2018_PDL1_treatment.Rda"
)

easierData <- data.frame(
  stringsAsFactors = FALSE,
  Title = "easierData",
  Description = paste(
    "Internal data required for proper performance of easier package. This data includes:",
    "cancer-specifc model parameters learned during training, mean and standard deviation of the",
    "expression of each gene across all TCGA cancer types, the list of genes used to compute",
    "correlated proxies of immune response, cancer-specific intercellular networks based on",
    "literature supported pairs from the Ramilowski database (Ramilowski et al., Nat.Commun., 2015),",
    "the frequency of each ligand-receptor pair feature across the whole TCGA database,",
    "information on how ligand-receptor pairs were grouped, annotations on approved gene symbols",
    "gene signatures included in the immune resistance program from Jerby-Arnon et al., 2018."
  ),
  BiocVersion = "3.13",
  Genome = NA,
  SourceType = "tar.gz",
  SourceUrl = NA,
  SourceVersion = NA,
  Species = "Homo sapiens",
  TaxonomyId = "9606",
  Coordinate_1_based = NA,
  DataProvider = NA,
  Maintainer = "Federico Marini <marinif@uni-mainz.de>",
  RDataClass = "list",
  DispatchClass = "Rda",
  RDataPath = "easier/easierData.Rda"
)

# RMTLR_TCGA_cancer_spec_trained_models <- data.frame(
#   stringsAsFactors = FALSE,
#   Title = "RMTLR_TCGA_cancer_spec_trained_models",
#   Description = paste(
#     "Cancer-specifc model parameters learned during training for different quantitative",
#     "descriptors using Regularized Multi-Task Linear Regression (RMTLR) algorithm with",
#     "randomized cross-validation repeated 100 times. A total of 100 models are available",
#     "for each task, here 10 tasks were used.",
#     "Also, information on training set statistics (i.e. meand and standard deviation)",
#     "are also provided for normalizing the test set accordingly.",
#     "Represented as a list object with both cancer-specific models and statistics about",
#     "the training dataset."
#   ),
#   BiocVersion = "3.13",
#   Genome = NA,
#   SourceType = "tar.gz",
#   SourceUrl = NA,
#   SourceVersion = NA,
#   Species = "Homo sapiens",
#   TaxonomyId = "9606",
#   Coordinate_1_based = NA,
#   DataProvider = NA,
#   Maintainer = "Federico Marini <marinif@uni-mainz.de>",
#   RDataClass = "list",
#   DispatchClass = "Rda",
#   RDataPath = "easier/RMTLR_TCGA_cancer_spec_trained_models.Rda"
# )
#
# TCGA_pancancer_gene_statistics <- data.frame(
#   stringsAsFactors = FALSE,
#   Title = "TCGA_pancancer_per_gene_statistics",
#   Description = paste(
#     "Mean and standard deviation of the expression of each gene across all TCGA cancer types.",
#     "Represented as a list object with both statistical measures. Helpful for normalization steps."
#   ),
#   BiocVersion = "3.13",
#   Genome = NA,
#   SourceType = "tar.gz",
#   SourceUrl = NA,
#   SourceVersion = NA,
#   Species = "Homo sapiens",
#   TaxonomyId = "9606",
#   Coordinate_1_based = NA,
#   DataProvider = NA,
#   Maintainer = "Federico Marini <marinif@uni-mainz.de>",
#   RDataClass = "list",
#   DispatchClass = "Rda",
#   RDataPath = "easier/TCGA_pancancer_gene_statistics.Rda"
# )
#
# proxies_genes <- data.frame(
#   stringsAsFactors = FALSE,
#   Title = "proxies_genes",
#   Description = paste(
#     "Genes used to compute correlated proxies of immune response."
#   ),
#   BiocVersion = "3.13",
#   Genome = NA,
#   SourceType = "tar.gz",
#   SourceUrl = NA,
#   SourceVersion = NA,
#   Species = "Homo sapiens",
#   TaxonomyId = "9606",
#   Coordinate_1_based = NA,
#   DataProvider = NA,
#   Maintainer = "Federico Marini <marinif@uni-mainz.de>",
#   RDataClass = "character",
#   DispatchClass = "Rda",
#   RDataPath = "easier/proxies_genes.Rda"
# )
#
# intercellular_cancer_spec_network <- data.frame(
#   stringsAsFactors = FALSE,
#   Title = "intercellular_cancer_spec_network",
#   Description = paste(
#     "Cancer-specific intercellular networks based on literature supported pairs",
#     "from the Ramilowski database (Ramilowski et al., Nat.Commun., 2015), filtering",
#     "for 24 cell types acknowledged to be present in the TME.",
#     "Additionally, a pan-cancer cell type is included using data from the CCLE",
#     "(Barretina et al., Nature, 2012) by computing the median expression of each",
#     "gene across cell lines related to the selected cancer types.",
#     "Represented as a list object with cancer-specific intercellular networks information."
#   ),
#   BiocVersion = "3.13",
#   Genome = NA,
#   SourceType = "tar.gz",
#   SourceUrl = NA,
#   SourceVersion = NA,
#   Species = "Homo sapiens",
#   TaxonomyId = "9606",
#   Coordinate_1_based = NA,
#   DataProvider = NA,
#   Maintainer = "Federico Marini <marinif@uni-mainz.de>",
#   RDataClass = "list",
#   DispatchClass = "Rda",
#   RDataPath = "easier/intercellular_cancer_spec_network.Rda"
# )
#
# LR_TCGA_frequency <- data.frame(
#   stringsAsFactors = FALSE,
#   Title = "LR_TCGA_frequency",
#   Description = paste(
#     "Frequency of each ligand-receptor pair across the whole TCGA database."
#   ),
#   BiocVersion = "3.13",
#   Genome = NA,
#   SourceType = "tar.gz",
#   SourceUrl = NA,
#   SourceVersion = NA,
#   Species = "Homo sapiens",
#   TaxonomyId = "9606",
#   Coordinate_1_based = NA,
#   DataProvider = NA,
#   Maintainer = "Federico Marini <marinif@uni-mainz.de>",
#   RDataClass = "numeric",
#   DispatchClass = "Rda",
#   RDataPath = "easier/LR_TCGA_frequency.Rda"
# )
#
# grouping_LR_pairs <- data.frame(
#   stringsAsFactors = FALSE,
#   Title = "grouping_LR_pairs",
#   Description = paste(
#     "Information on how ligand-receptor pairs were grouped because of sharing",
#     "the same gene either as ligand or receptor.",
#     "Represented as a list object with information to group certain LR pairs."
#   ),
#   BiocVersion = "3.13",
#   Genome = NA,
#   SourceType = "tar.gz",
#   SourceUrl = NA,
#   SourceVersion = NA,
#   Species = "Homo sapiens",
#   TaxonomyId = "9606",
#   Coordinate_1_based = NA,
#   DataProvider = NA,
#   Maintainer = "Federico Marini <marinif@uni-mainz.de>",
#   RDataClass = "list",
#   DispatchClass = "Rda",
#   RDataPath = "easier/grouping_LR_pairs.Rda"
# )
#
# HGNC_checker <- data.frame(
#   stringsAsFactors = FALSE,
#   Title = "HGNC_checker",
#   Description = paste(
#     "Information on gene symbols approved annotations.",
#     "Represented as a matrix, derived from https://www.genenames.org/tools/multi-symbol-checker/."
#   ),
#   BiocVersion = "3.13",
#   Genome = NA,
#   SourceType = "tar.gz",
#   SourceUrl = "https://www.genenames.org/tools/multi-symbol-checker",
#   SourceVersion = NA,
#   Species = "Homo sapiens",
#   TaxonomyId = "9606",
#   Coordinate_1_based = NA,
#   DataProvider = NA,
#   Maintainer = "Federico Marini <marinif@uni-mainz.de>",
#   RDataClass = "data.frame",
#   DispatchClass = "Rda",
#   RDataPath = "easier/HGNC_checker.Rda"
# )
#
# RIR_signature <- data.frame(
#   stringsAsFactors = FALSE,
#   Title = "RIR_signature",
#   Description = paste(
#     "Gene signatures included in the immune resistance program from Jerby-Arnon et al., 2018",
#     "Represented as a list object obtained from",
#     "https://github.com/livnatje/ImmuneResistance/blob/master/Results/Signatures/resistance.program.RData."
#   ),
#   BiocVersion = "3.13",
#   Genome = NA,
#   SourceType = "tar.gz",
#   SourceUrl = "https://github.com/livnatje/ImmuneResistance/blob/master/Results/Signatures/resistance.program.RData",
#   SourceVersion = NA,
#   Species = "Homo sapiens",
#   TaxonomyId = "9606",
#   Coordinate_1_based = NA,
#   DataProvider = NA,
#   Maintainer = "Federico Marini <marinif@uni-mainz.de>",
#   RDataClass = "list",
#   DispatchClass = "Rda",
#   RDataPath = "easier/RIR_signature.Rda"
# )


# write to .csv
write.csv(rbind(
  Mariathasan2018_PDL1_treatment,
  easierData
), file = "inst/extdata/metadata.csv", row.names = FALSE)
