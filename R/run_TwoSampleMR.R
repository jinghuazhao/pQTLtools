#' Basic TwoSampleMR analysis
#'
#' Given harmonised data, this function conducts a two-sample MR analysis.
#'
#' @param TwoSampleMRinput Harmonised data.
#' @param mr_plot one of "None", "TwoSampleMR", "pQTLtools" for no, the original and the revised plots, respectively.
#' @param prefix a prefix for output files.
#'
#' @details
#' As TwoSampleMR faces seemingly perplexing options, this function intends to simplify various steps in a two-sample MR
#' as in \insertCite{dt18;textual}{pQTLtools}. It is particularly useful when a large numbher of MRs are necessary,
#' e.g., multiple proteins and their cis/trans regions need to be examined, in which case prefix could direct the output
#' to various directories.
#'
#' @export
#' @return
#' No value is returned but several files.
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' library(TwoSampleMR)
#' library(pQTLtools)
#' prot <- "MMP.10"
#' type <- "cis"
#' f <- paste0(prot,"-",type,".mrx")
#' d <- read.table(file.path(find.package("pQTLtools",lib.loc=.libPaths()),"tests",f),
#'                 header=TRUE)
#' exposure <- TwoSampleMR::format_data(within(d,{P=10^logP}), phenotype_col="prot", snp_col="rsid",
#'                                      chr_col="Chromosome", pos_col="Posistion",
#'                                      effect_allele_col="Allele1", other_allele_col="Allele2",
#'                                      eaf_col="Freq1", beta_col="Effect", se_col="StdErr",
#'                                      pval_col="P", log_pval=FALSE,
#'                                      samplesize_col="N")
#' clump <- exposure[sample(1:nrow(exposure),nrow(exposure)/80),] # TwoSampleMR::clump_data(exposure)
#' outcome <- TwoSampleMR::extract_outcome_data(snps=exposure$SNP,outcomes="ebi-a-GCST007432")
#' harmonise <- TwoSampleMR::harmonise_data(clump,outcome)
#' prefix <- paste(prot,type,sep="-")
#' run_TwoSampleMR(harmonise, mr_plot="pQTLtools", prefix=prefix)
#' caption <- "Table. MMP.10 variants and FEV1"
#' knitr::kable(read.delim(paste0(prefix,"-result.txt"),header=TRUE),
#'              caption=paste(caption, "(result)"))
#' knitr::kable(read.delim(paste0(prefix,"-heterogeneity.txt"),header=TRUE),
#'              caption=paste(caption,"(heterogeneity)"))
#' knitr::kable(read.delim(paste0(prefix,"-pleiotropy.txt"),header=TRUE),
#'              caption=paste(caption,"(pleiotropy)"))
#' knitr::kable(read.delim(paste0(prefix,"-single.txt"),header=TRUE),
#'              caption=paste(caption,"(single)"))
#' knitr::kable(read.delim(paste0(prefix,"-loo.txt"),header=TRUE),
#'              caption=paste(caption,"(loo)"))
#' for (x in c("result","heterogeneity","pleiotropy","single","loo"))
#'     unlink(paste0(prefix,"-",x,".txt"))
#'
#' @keywords utilities

run_TwoSampleMR <- function(TwoSampleMRinput, mr_plot="None", prefix="")
{
  for (p in "TwoSampleMR") {
    if (length(grep(paste("^package:", p, "$", sep=""), search())) == 0) {
       if (!requireNamespace(p, quietly = TRUE))
           warning(paste("This function needs package `", p, "' to be fully functional; please install", sep=""))
    }
  }
  harmonise <- TwoSampleMRinput
  result <- heterogeneity <- pleiotropy <- single <- loo <- NULL
  result <- TwoSampleMR::mr(harmonise)
  heterogeneity <- TwoSampleMR::mr_heterogeneity(harmonise)
  pleiotropy <- TwoSampleMR::mr_pleiotropy_test(harmonise)
  single <- TwoSampleMR::mr_singlesnp(harmonise)
  loo <- TwoSampleMR::mr_leaveoneout(harmonise)
  invisible(lapply(c("result","heterogeneity","pleiotropy","single","loo"), function(x)
                   if (exists(x)) write.table(format(get(x),digits=3),
                                              file=paste0(prefix,"-",x,".txt"),
                                              quote=FALSE,row.names=FALSE,sep="\t")))
  type <- match.arg(mr_plot, c("None", "TwoSampleMR", "pQTLtools"))
  if (type == 1) return(0) else if (type == "TwoSampleMR")
  { 
    scatter = TwoSampleMR::mr_scatter_plot(result, harmonise)
    forest <- TwoSampleMR::mr_forest_plot(single)
    funnel <- TwoSampleMR::mr_funnel_plot(single)
    leaveoneout <- TwoSampleMR::mr_leaveoneout_plot(loo)
  } else {
    scatter = mr_scatter_plot2(result, harmonise)
    forest <- mr_forest_plot2(single)
    funnel <- mr_funnel_plot2(single)
    leaveoneout <- mr_leaveoneout_plot2(loo)
  }
  invisible(sapply(c(scatter,forest,funnel,leaveoneout), function(x) print(x)))
}
