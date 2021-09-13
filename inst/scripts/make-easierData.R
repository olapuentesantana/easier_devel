# Organize internal data needed by easier for proper performance

# *******************************************
##  RMTLR_TCGA_cancer_spec_trained_models ##
# *******************************************
# Cancer-specifc model parameters learned during training for different quantitative
# descriptors using Regularized Multi-Task Linear Regression (RMTLR) algorithm with
# randomized cross-validation repeated 100 times. A total of 100 models are available
# for each task, here 10 tasks were used.
# Also, information on training set statistics (i.e. meand and standard deviation
# are also provided for normalizing the test set accordingly.
# Models and training set statistics represented as a list.

opt_models
opt_xtrain_stats


# *******************************************
##  TCGA_pancancer_gene_statistics ##
# *******************************************

# Mean and standard deviation of the expression of each gene across all TCGA cancer types.
# Each measure represented as a numeric vector Helpful for normalization steps.

TCGA_mean_pancancer
TCGA_sd_pancancer


# *******************************************
##  proxies_genes ##
# *******************************************

# Genes used to compute correlated proxies of immune response

cor_genes_to_remove

# *******************************************
##  intercellular_cancer_spec_network ##
# *******************************************

# Cancer-specific intercellular networks based on literature supported pairs
# from the Ramilowski database (Ramilowski et al., Nat.Commun., 2015), filtering
# for 24 cell types acknowledged to be present in the TME.
# Additionally, a pan-cancer cell type is included using data from the CCLE
# (Barretina et al., Nature, 2012) by computing the median expression of each
# gene across cell lines related to the selected cancer types.
# Represented as a list object with cancer-specific intercellular networks information.

intercell_network_cancer_spec

# *******************************************
##  LR_TCGA_frequency ##
# *******************************************

# Frequency of each ligand-receptor pair across the whole TCGA database.

lr_frequency

# *******************************************
##  grouping_LR_pairs ##
# *******************************************

# Information on how ligand-receptor pairs were grouped because of sharing
# the same gene either as ligand or receptor.
# Represented as a list object with information to group certain LR pairs.

grouping_lrpairs_info

# *******************************************
##  HGNC_checker ##
# *******************************************

# Information on gene symbols approved annotations.
# Represented as a matrix, derived from https://www.genenames.org/tools/multi-symbol-checker/.

HGNC

# *******************************************
##  RIR_signature ##
# *******************************************

# Gene signatures included in the immune resistance program from Jerby-Arnon et al., 2018.
# Represented as a list object obtained from
# https://github.com/livnatje/ImmuneResistance/blob/master/Results/Signatures/resistance.program.RData."

res_sig

# save easier data
save(opt_models,
  opt_xtrain_stats,
  TCGA_mean_pancancer,
  TCGA_sd_pancancer,
  cor_genes_to_remove,
  intercell_network_cancer_spec,
  lr_frequency,
  grouping_lrpairs_info,
  HGNC,
  res_sig,
  file = "inst/extdata/easierData.rda"
)
