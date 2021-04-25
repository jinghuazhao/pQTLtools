pqtlMR <- function(Ins,ids,outdir="./")
{
# 0. library(MRInstruments)
  library(TwoSampleMR)
# 1. read in the exposure data from a file
  Ins <- format_data(Ins, type = "exposure", header = TRUE,
                   phenotype_col = "Phenotype", snp_col = "SNP", beta_col = "beta",
                   se_col = "se", eaf_col = "eaf", effect_allele_col = "effect_allele",
                   other_allele_col = "other_allele", pval_col = "pval")
# 2. read in the outcome data
# ao <- available_outcomes(access_token=NULL)
  outcome_dat <- extract_outcome_data(snps = Ins$SNP,outcomes = ids)
# 3. harmonise the exposure and outcome data
  dat <- NULL
  dat <- harmonise_data(exposure_dat = Ins, outcome_dat = outcome_dat)
# 4. run the MR and sensitivity analyses 
  mr_results <- NULL # main MR analysis
  mr_hetero <- NULL  # heterogeneity test across instruments
  mr_pleio <- NULL   # MR-Egger intercept test
  mr_single <- NULL  # single SNP MR using Wald ratio
  try(mr_results <- mr(dat, method_list=c("mr_wald_ratio", "mr_ivw")))
  mr_hetero <- mr_heterogeneity(dat)
  mr_pleio <- mr_pleiotropy_test(dat)
  try(mr_single <- mr_singlesnp(dat))
# 5. save the MR results
  harmonise <- dat
  invisible(lapply(c("harmonise","mr_results","mr_hetero","mr_pleio","mr_single"),
                   function(x) if(exists(x)) write.table(format(get(x),digits=3),
                               file=file.path(outdir,paste0("pQTL-combined-",x,".txt")),
                               col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")))
}
