#===============================================================================
# twosample MR analysis using MR-base
#===============================================================================

rm(list=ls())

# install.packages("remotes")
# remotes::install_github("MRCIEU/TwoSampleMR")
library(TwoSampleMR)

setwd("G:/UKBiobank_AnnualReport_urgent!/Iqbal/madakkatel2024large/data/interim")

## load the MRbase outcome database
ao <- available_outcomes()
colnames(ao)

# Vector list of exposure selected (included) in the MR
list_of_studies <- c(
  "ukb-b-17422", "ieu-a-1004", "ukb-b-1209", "ieu-b-4760", "ebi-a-GCST90000050",
  "finn-b-N14_FEMALEINFERT", "met-d-Omega_3", "met-c-855", "met-d-Omega_6",
  "met-c-856", "met-d-Omega_6_by_Omega_3", "ukb-b-3460", "ieu-b-73",
  "ukb-b-16998", "ebi-a-GCST006099", "ebi-a-GCST006097", "ebi-a-GCST006098",
  "ieu-a-1187", "ukb-b-19809", "ebi-a-GCST006947", "ukb-a-35", "ukb-b-10787",
  "ieu-a-89", "ukb-b-16881", "ieu-a-1070", "ukb-b-4650", "ieu-a-1096",
  "ukb-b-11842", "ieu-a-107", "ukb-b-19393", "ieu-a-999", "ukb-b-20044",
  "ukb-b-6704", "ukb-b-8338", "ukb-b-18096", "ukb-b-7212", "ukb-b-16446",
  "ukb-b-8909", "ebi-a-GCST003435", "ukb-b-16407", "ukb-b-12854", "ukb-b-20188",
  "ukb-b-20531", "ukb-b-18377", "ukb-b-9405", "ieu-a-61", "ukb-b-15590",
  "ieu-a-49", "ukb-b-20175", "ieu-b-4818", "ukb-b-7992", "ieu-b-39",
  "ukb-b-15892", "ukb-b-10215", "ukb-b-7478", "ukb-b-3376", "ukb-b-12019",
  "ukb-b-7953", "ieu-b-105", "ukb-b-19657", "ebi-a-GCST007432",
  "ebi-a-GCST007431", "ukb-d-30650_irnt", "ukb-d-30620_irnt", "ukb-d-30770_irnt",
  "prot-a-1443", "ukb-d-30710_irnt", "ieu-b-4764", "ukb-d-30690_irnt",
  "met-a-307", "ieu-b-109", "ebi-a-GCST002223", "ieu-b-111", "met-c-934",
  "ukb-d-30740_irnt", "ebi-a-GCST90002232", "ukb-d-30750_irnt",
  "ukb-d-30720_irnt", "prot-c-2609_59_2", "ukb-d-30860_irnt", "ukb-d-30070_irnt",
  "ebi-a-GCST006804", "ebi-a-GCST90002404", "ukb-d-30040_irnt",
  "ebi-a-GCST004602", "ukb-d-30050_irnt", "ukb-d-30290_irnt", "ukb-d-30270_irnt",
  "ukb-d-30120_irnt", "ukb-d-30180_irnt", "ukb-d-30140_irnt", "ieu-b-34",
  "ukb-d-30200_irnt", "ebi-a-GCST004606", "ukb-d-30210_irnt", "ukb-d-30190_irnt",
  "ukb-d-30130_irnt", "ieu-b-31"
)
# extract variant-exposure data
exposure_selected <- ao[ao$id %in%list_of_studies,]
rownames(exposure_selected) <- NULL
variant_exposure <- extract_instruments(outcomes=exposure_selected$id,
                                        p1 = 5e-08,
                                        clump = TRUE,
                                        p2 = 5e-08,
                                        r2 = 0.001,
                                        kb = 10000,
                                        access_token = ieugwasr::check_access_token(),
                                        force_server = FALSE
)

# extract variant-outcome data
idoutcome <- c('ieu-a-1120','ieu-a-1228','ieu-a-1125', 'ieu-a-1124','ieu-a-1123')
outcome_selected <- ao[ao$id%in%idoutcome,]
rownames(outcome_selected) <- NULL

variant_outcome <- extract_outcome_data(snps = variant_exposure$SNP,
                                        outcomes = idoutcome,
                                        proxies = TRUE,
                                        rsq = 0.8,
                                        align_alleles = 1,
                                        palindromes = 1,
                                        maf_threshold = 0.3)

# harmonization of exposure and outcome data
dat <- harmonise_data(variant_exposure, variant_outcome, action = 2)

# MR-analysis
mr_results <- TwoSampleMR::mr(dat, method_list=c("mr_wald_ratio","mr_ivw","mr_egger_regression"))

# Restrict the R environment
keep <- c("list_of_studies","dat","mr_results", "exposure_selected", "outcome_selected")
to_remove <- setdiff(ls(), keep)
rm(list = to_remove)

# save the R environment
save.image("G:/UKBiobank_AnnualReport_urgent!/Iqbal/madakkatel2024large/data/interim/mr_results_mr-base_all.RData")

