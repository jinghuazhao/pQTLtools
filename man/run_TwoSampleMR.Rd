\name{run_TwoSampleMR}
\alias{run_TwoSampleMR}
\title{A generic wrapper for TwoSampleMR analysis}
\usage{run_TwoSampleMR(TwoSampleMRinput, plot=TRUE, prefix)}
\description{
Given a compressed or uncompressed file containing exposure data, this function conducts a two-sample MR analysis using MR-Base.
}
\arguments{
\item{TwoSampleMRinput}{A list containing the following elements to allow for greater flexibility,
  \enumerate{
   \item{exposure}{result TwoSampleMR's format_data()}
   \item{outcome}{result if TwoSampleMR's extract_outcome_data()}
   \item{clump}{result of TwoSampleMR's clump_data()}
   \item{harmonise}{result of TwoSampleMR's harmonise_data()}
  }
}
\item{plot}{to produce plots}
\item{prefix}{a prefix for output files}
}
\details{
As TwoSampleMR includes seemingly daunting tasks, this function intends to simplify various steps in a two-sample MR. It is
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
outcomes <- "ebi-a-GCST007432"
prot <- "MMP.10"
type <- "cis"
f <- paste0(prot,"-",type,".mrx")
d <- read.table(file.path(path.package("pQTLtools"),"tests",f),header=TRUE)
exposure <- format_data(within(d,{P=10^logP}), phenotype_col="prot", snp_col="rsid",
                        chr_col="Chromosome", pos_col="Posistion",
                        effect_allele_col="Allele1", other_allele_col="Allele2",
                        eaf_col="Freq1", beta_col="Effect", se_col="StdErr",
                        pval_col="P", log_pval=FALSE,
                        samplesize_col="N")
clump <- clump_data(exposure)
outcome <- extract_outcome_data(snps=exposure$SNP,outcomes=outcomes)
harmonise <- harmonise_data(clump,outcome)
input <- list(exposure,outcome,clump,harmonise)
run_TwoSampleMR(input, prefix=paste(outcomes,prot,type,sep="-"))
}
\note{Adapted from script by Dimou NL, Tsilidis KK}
\keyword{utilities}
