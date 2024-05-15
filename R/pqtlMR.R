#' Basic pQTL-MR analysis
#'
#' This function takes exposure and outcome data as produced by `format_data()` and used to perform
#' MR analysis against a list of outcomes; in the latter case it can be data from MR-Base, e.g.,
#' `outcome <- extract_outcome_data(snps=with(exposure,SNP),outcomes=c("ieu-a-7","ebi-a-GCST007432"))`.
#'
#' @param exposure exposure data.
#' @param outcome the counterpart for outcome.
#' @param mr_plot to produce plots.
#' @param prefix a prefix for output files.
#' @param reverse if TRUE, perform reverse MR.
#'
#' @details
#' Adapted from \insertCite{zheng20;textual}{pQTLtools}, this function is analogous to `run_TwoSampleMR()`.
#' @md
#'
#' @export
#' @return
#' No value is returned but several files.
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' fi <- file.path(find.package("pQTLtools",lib.loc=.libPaths()),"tests","Ins.csv")
#' exposure <- TwoSampleMR::format_data(read.csv(fi),type="exposure")
#' fo <- file.path(find.package("pQTLtools",lib.loc=.libPaths()),"tests","Out.csv")
#' outcome <- TwoSampleMR::format_data(read.csv(fo),type="outcome")
#' pQTLtools::pqtlMR(exposure, outcome, prefix="IL6R-")
#' pQTLtools::pqtlMR(exposure, outcome, prefix="IL6R_rev-",reverse=TRUE)
#' unlink(c("IL6R*","pQTL-combined*"))
#' # Phenotype,SNP,effect_allele,other_allele,eaf,beta,se,pval
#' # ABO,rs505922,C,T,0.313,1.298,0.014,1.2e-1828
#' # LIFR,rs635634,T,C,0.180,-0.300,0.032,6.00E-21
#' # f <- file.path(find.package("pQTLtools",lib.loc=.libPaths()),"tests","ms.ins")
#' # exposure <- format_data(read.table(f, header=TRUE), samplesize_col="N")
#' # SNP Phenotype effect_allele other_allele eaf beta se pval N
#' # rs1800693 TNFB T C 0.6033 0.0282  0.0136 0.0389045   11787
#' # rs2364485 TNFB A C 0.1645 6514963 0.1759 1.62181e-20 11344
#' # https://raw.githubusercontent.com/MRCIEU/epigraphdb-pqtl/master/scripts/MR-pQTL-script.R
#'
#' @note
#' Adapted from script by Jie Zheng.
#' @keywords utilities

pqtlMR <- function(exposure, outcome, mr_plot=FALSE, prefix="pQTL-combined-", reverse=FALSE)
{
   for (p in c("TwoSampleMR")) {
     if (length(grep(paste("^package:", p, "$", sep=""), search())) == 0) {
        if (!requireNamespace(p, quietly = TRUE))
            warning(paste("This function needs package `", p, "' to be fully functional; please install", sep=""))
     }
   }
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
   harmonise <- TwoSampleMR::harmonise_data(exposure,outcome)
   if (reverse) harmonise <- subset(within(harmonise,
   {
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
   result <- single <- NULL
   try(result <- TwoSampleMR::mr(harmonise,method_list=c("mr_wald_ratio", "mr_ivw")))
   try(single <- TwoSampleMR::mr_singlesnp(harmonise))
   invisible(lapply(c("harmonise","result","single"), function(x) if (exists(x))
             write.table(format(get(x),digits=3),file=paste0(prefix,x,".txt"),
                         quote=FALSE,row.names=FALSE,sep="\t")))
   if (mr_plot)
   {
      scatter <- TwoSampleMR::mr_scatter_plot(result, harmonise)
      forest <- TwoSampleMR::mr_forest_plot(single)
      funnel <- TwoSampleMR::mr_funnel_plot(single)
      invisible(sapply(c(scatter,forest,funnel), function(x) print(x)))
   }
}
