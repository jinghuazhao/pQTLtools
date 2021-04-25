format_file.args <- function(file=NULL, type = "exposure", snps = NULL, header = TRUE,
                             phenotype_col = "Phenotype", snp_col = "SNP", beta_col = "beta",
                             se_col = "se", eaf_col = "eaf", effect_allele_col = "effect_allele",
                             other_allele_col = "other_allele", pval_col = "pval", units_col = "units",
                             ncase_col = "ncase", ncontrol_col = "ncontrol", samplesize_col = "samplesize",
                             gene_col = "gene", id_col = "id", min_pval = 1e-200, z_col = "z",
                             info_col = "info", chr_col = "chr", pos_col = "pos", log_pval = FALSE)
  invisible(list(file = file, type = type, snps = snps, header = header,
                 phenotype_col = phenotype_col, snp_col = snp_col, beta_col = beta_col,
                 se_col = se_col, eaf_col = eaf_col, effect_allele_col = effect_allele_col,
                 other_allele_col = other_allele_col, pval_col = pval_col, units_col = units_col,
                 ncase_col = ncase_col, ncontrol_col = ncontrol_col, samplesize_col = samplesize_col,
                 gene_col = gene_col, id_col = id_col, min_pval = min_pval, z_col = z_col,
                 info_col = info_col, chr_col = chr_col, pos_col = pos_col, log_pval = FALSE))

extract_outcome_data.args <- function(snps=NULL, outcomes=NULL, proxies=TRUE, rsq=0.8, align_alleles=1,
                                      palindromes=1, maf_threshold=0.3, access_token=ieugwasr::check_access_token(),
                                      splitsize=10000, proxy_splitsize=500)
  invisible(list(snps=snps,outcomes=outcomes,proxies=proxies,rsq=rsq,align_alleles=align_alleles,
                           palindromes=palindromes,maf_threshold=maf_threshold,access_token=access_token,
                           splitsize=splitsize,proxy_splitsize=proxy_splitsize))

clump_data.args <- function(dat=NULL, clump_kb = 10000, clump_r2 = 0.001, clump_p1 = 1, clump_p2 = 1, pop = "EUR")
  invisible(list(dat=dat,clump_kb=clump_kb,clump_r2=clump_r2,clump_p1=clump_p1,clump_p2=clump_p2,pop=pop))

harmonise_data.args <- function(exposure_dat=NULL, outcome_dat=NULL, action = 2)
  invisible(list(exporesure_dat=exposure_dat,outcome_dat=outcome_dat,action=action))

swap <- function(x,y)
   eval(parse(text = paste("swap_unique_var_a <-", substitute(x), ";",
   substitute(x), "<-", substitute(y), ";",
   substitute(y), "<-swap_unique_var_a")), envir=parent.frame())

pqtlMR <- function(Ins=format_file.args(),Ids=extract_outcome_data.args(),harmonise=harmonise_data.args(),
                   prefix="INF1",reverse=FALSE,...)
{
  exposure <- NA
  outcome <- NA
  id.exposure <- NA
  id.outcome <- NA
  effect_allele.exposure <- NA
  effect_allele.outcome <- NA
  other_allele.exposure <- NA
  other_allele.outcome <- NA
  eaf.exposure <- NA
  eaf.outcome <- NA
  samplesize.exposure <- NA
  beta.exposure <- NA
  beta.outcome <- NA
  se.exposure <- NA
  se.outcome <- NA
  pval.exposure <- NA
  pval.outcome <- NA
  swap_unique_var_a <- NA
  d <- with(Ins,lapply(file, function(x) tryCatch(read.table(file,as.is=TRUE,header=TRUE), error=function(e) NULL))[[1]])
  if (nrow(d)==0) stop("the Instrument data is empty")
  exposure_dat <- with(Ins,TwoSampleMR::format_data(d, type = type, snps = snps, header = header,
                       phenotype_col = phenotype_col, snp_col = snp_col, beta_col = beta_col,
                       se_col = se_col, eaf_col = eaf_col, effect_allele_col = effect_allele_col,
                       other_allele_col = other_allele_col, pval_col = pval_col, units_col = units_col,
                       ncase_col = ncase_col, ncontrol_col = ncontrol_col, samplesize_col = samplesize_col,
                       gene_col = gene_col, id_col = id_col, min_pval = min_pval, z_col = z_col,
                       info_col = info_col, chr_col = chr_col, pos_col = pos_col, log_pval = FALSE))
# ao <- TwoSampleMR::available_outcomes(access_token=NULL)
  if (is.null(Ids$snps)) Ids$snps <- exposure_dat$SNP
  outcome_dat <- with(Ids,TwoSampleMR::extract_outcome_data(snps, outcomes, proxies=proxies, rsq=rsq,
                      align_alleles=align_alleles, palindromes=palindromes, maf_threshold=maf_threshold))
  if (is.null(harmonise$exposure_dat)) harmonise$exposure_dat <- exposure_dat
  if (is.null(harmonise$outcome_dat)) harmonise$outcome_dat <- outcome_dat
  harmonise <- with(harmonise,TwoSampleMR::harmonise_data(exposure_dat, outcome_dat, action=action))
  if (reverse) harmonise <- subset(within(harmonise,
  {
    swap(exposure,outcome)
    swap(id.exposure,id.outcome)
    swap(effect_allele.exposure,effect_allele.outcome)
    swap(other_allele.exposure,other_allele.outcome)
    swap(eaf.exposure,eaf.outcome)
    if(exists("samplesize.exposure") & exists("samplesize.outcome")) swap(samplesize.exposure,samplesize.outcome)
    if(!exists("samplesize.exposure") & exists("samplesize.outcome")) samplesize.outcome <- NA
    swap(beta.exposure,beta.outcome)
    swap(se.exposure,se.outcome)
    swap(pval.exposure,pval.outcome)
  }),select=-swap_unique_var_a)
  result <- heterogeneity <- pleiotropy <- single <- NULL
  try(result <- TwoSampleMR::mr(harmonise, method_list=c("mr_wald_ratio", "mr_ivw"))) # main MR analysis
  heterogeneity <- TwoSampleMR::mr_heterogeneity(harmonise) # heterogeneity test across instruments
  pleiotropy <- TwoSampleMR::mr_pleiotropy_test(harmonise) # MR-Egger intercept test
  try(single <- TwoSampleMR::mr_singlesnp(harmonise)) #single SNP MR using Wald ratio
# result <- within(result,outcome <- sub(" [|]* id:ieu-a-[0-9]*| [|]* id:ukb-a-[0-9]*", "\\1", outcome, perl = TRUE))
  ext <- ".txt"
  invisible(lapply(c("harmonise","result","heterogeneity","pleiotropy","single"), function(x) {
                   v <- lapply(x, function(x) tryCatch(get(x), error=function(e) NULL))[[1]]
                   if (!is.null(v)) write.table(format(v,digits=3),file=paste0(prefix,"-",x,ext),quote=FALSE,row.names=FALSE,sep="\t")
            })
  )
}

run_TwoSampleMR <- function(exposure.args=format_file.args(),outcome.args=extract_outcome_data.args(),
                            clump.args=clump_data.args(),harmonise.args=harmonise_data.args(),prefix,...)
{
  d <- with(exposure.args,lapply(file, function(x) tryCatch(read.table(file,as.is=TRUE,header=TRUE), error=function(e) NULL))[[1]])
  if (nrow(d)==0) stop("the exposure data is empty")
  e <- with(exposure.args,TwoSampleMR::format_data(d, type = type, snps = snps, header = header,
                 phenotype_col = phenotype_col, snp_col = snp_col, beta_col = beta_col,
                 se_col = se_col, eaf_col = eaf_col, effect_allele_col = effect_allele_col,
                 other_allele_col = other_allele_col, pval_col = pval_col, units_col = units_col,
                 ncase_col = ncase_col, ncontrol_col = ncontrol_col, samplesize_col = samplesize_col,
                 gene_col = gene_col, id_col = id_col, min_pval = min_pval, z_col = z_col,
                 info_col = info_col, chr_col = chr_col, pos_col = pos_col, log_pval = FALSE))
  if (is.null(clump.args$dat)) clump.args$dat <- e
  exposure_dat <- with(clump.args,TwoSampleMR::clump_data(dat,clump_kb=clump_kb, clump_r2=clump_r2, clump_p1=clump_p1, clump_p2=clump_p2, pop=pop))
  if (is.null(outcome.args$snps)) outcome.args$snps <- with(e,SNP)
  outcome_dat <- with(outcome.args,TwoSampleMR::extract_outcome_data(snps, outcomes, proxies=proxies, rsq=rsq,
                                   align_alleles=align_alleles, palindromes=palindromes, maf_threshold=maf_threshold))
  if(!is.null(outcome_dat))
  {
    if (is.null(harmonise.args$exposure_dat)) harmonise.args$exposure_dat <- exposure_dat
    if (is.null(harmonise.args$outcome_dat)) harmonise.args$outcome_dat <- outcome_dat
    dat <- with(harmonise.args,TwoSampleMR::harmonise_data(exposure_dat, outcome_dat, action=action))
    TwoSampleMR::directionality_test(dat)
    if (nrow(dat)!=0)
    {
      result <- TwoSampleMR::mr(dat)
      heterogeneity <- TwoSampleMR::mr_heterogeneity(dat)
      pleiotropy <- TwoSampleMR::mr_pleiotropy_test(dat)
      single <- TwoSampleMR::mr_singlesnp(dat)
      loo <- TwoSampleMR::mr_leaveoneout(dat)
      invisible(lapply(c("result","heterogeneity","pleiotropy","single","loo"), function(x) {
                      v <- lapply(x, function(x) tryCatch(get(x), error=function(e) NULL))[[1]]
                      if (!is.null(v)) write.table(format(v,digits=3),file=paste0(prefix,"-",x,".txt"),quote=FALSE,row.names=FALSE,sep="\t")
                    }))
      scatter <- TwoSampleMR::mr_scatter_plot(result, dat)
      forest <- TwoSampleMR::mr_forest_plot(single)
      funnel <- TwoSampleMR::mr_funnel_plot(single)
      leaveoneout <- TwoSampleMR::mr_leaveoneout_plot(loo)
      pdf(paste0(prefix,".pdf"))
      invisible(sapply(c(scatter,forest,funnel,leaveoneout), function(x) print(x)))
      dev.off()
    }
  }
}
