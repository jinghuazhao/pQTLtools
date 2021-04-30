\name{run_TwoSampleMR}
\alias{run_TwoSampleMR}
\title{Basic TwoSampleMR analysis}
\usage{run_TwoSampleMR(TwoSampleMRinput, mr_plot="None", prefix)}
\description{
Given harmonised data, this function conducts a two-sample MR analysis.
}
\arguments{
\item{TwoSampleMRinput}{Harmonised data}
\item{mr_plot}{one of "None", "TwoSampleMR", "New" for no, the original and (some) new plots, respectively}
\item{prefix}{a prefix for output files}
}
\details{
As TwoSampleMR faces seemingly perplexing options, this function intends to simplify various steps in a two-sample MR. It is
particularly useful when a large numbher of MRs are necessary, e.g., multiple proteins and their cis/trans regions need to be examined,
in which case prefix could direct the output to various directories.
}

\value{
No value is returned but several files.
}

\section{References}{
Dimou NL, Tsilidis KK. A Primer in Mendelian Randomization Methodology with a Focus on Utilizing Published Summary Association Data. In
Evangelos Evangelou (ed.), Genetic Epidemiology: Methods and Protocols, Methods in Molecular Biology, vol. 1793,
https://doi.org/10.1007/978-1-4939-7868-7_13, Springer Science+Business Media, LLC, part of Springer Nature 2018
}

\examples{
library(TwoSampleMR)
library(pQTLtools)
outcomes <- "ebi-a-GCST007432"
prot <- "MMP.10"
type <- "cis"
f <- paste0(prot,"-",type,".mrx")
d <- read.table(file.path(find.package("pQTLtools",lib.loc=.libPaths()),"tests",f),
                header=TRUE)
exposure <- format_data(within(d,{P=10^logP}), phenotype_col="prot", snp_col="rsid",
                        chr_col="Chromosome", pos_col="Posistion",
                        effect_allele_col="Allele1", other_allele_col="Allele2",
                        eaf_col="Freq1", beta_col="Effect", se_col="StdErr",
                        pval_col="P", log_pval=FALSE,
                        samplesize_col="N")
clump <- clump_data(exposure)
outcome <- extract_outcome_data(snps=exposure$SNP,outcomes=outcomes)
harmonise <- harmonise_data(clump,outcome)
prefix <- paste(outcomes,prot,type,sep="-")
run_TwoSampleMR(harmonise, mr_plot="New", prefix=prefix)
caption <- "Table. MMP.10 variants and FEV1"
knitr::kable(read.delim(paste0(prefix,"-result.txt"),header=TRUE),
             caption=paste(caption, "(result)"))
knitr::kable(read.delim(paste0(prefix,"-heterogeneity.txt"),header=TRUE),
             caption=paste(caption,"(heterogeneity)"))
knitr::kable(read.delim(paste0(prefix,"-pleiotropy.txt"),header=TRUE),
             caption=paste(caption,"(pleiotropy)"))
knitr::kable(read.delim(paste0(prefix,"-single.txt"),header=TRUE),
             caption=paste(caption,"(single)"))
knitr::kable(read.delim(paste0(prefix,"-loo.txt"),header=TRUE),
             caption=paste(caption,"(loo)"))
for (x in c("result","heterogeneity","pleiotropy","single","loo"))
    unlink(paste0(prefix,"-",x,".txt"))
}
\note{Adapted from script by Dimou NL, Tsilidis KK}
\keyword{utilities}
