\name{run_TwoSampleMR}
\alias{run_TwoSampleMR}
\title{A generic wrapper for TwoSampleMR analysis}
\usage{run_TwoSampleMR(exposure.args=format_data.args(),outcome.args=extract_outcome_data.args(),
                       clump.args=clump_data.args(),harmonise.args=harmonise_data.args(),prefix,...)}
\description{
Given a compressed or uncompressed file containing exposure data, this function conducts a two-sample MR analysis using MR-Base.
}
\arguments{
\item{exposure.args}{a list of options appropriate for TwoSampleMR's format_data(),
  \enumerate{
    \item file. file
    \item type. exposure or outcome
    \item phenotype_col. protein column
    \item header. TRUE or FALSE
    \item snp_col. rsid
    \item effect_allele_col. effect allele
    \item other_allele_col. other allele
    \item eaf_col. EAF
    \item beta_col. beta
    \item se_col. se
    \item pval_col. pval
    \item log_pval. FALSE or TRUE
  }
}
\item{outcome.args}{a list of options appropriate for TwoSampleMR's extract_outcome_data(),
  \enumerate{
    \item snp. NULL (to derive from exposure.args)
    \item outcomes. NULL
    \item proxies. TRUE or FALSE
    \item rsq. r^2
    \item align_alleles. 1
    \item palindromes. 1
    \item maf_threshold. 0.3
    \item access_token. token
    \item splitsize. 10000
    \item proxy_splitsize. 500
  }
}
\item{clump.args}{a list of options appropriate for TwoSampleMR's clump_data(),
  \enumerate{
    \item dat. NULL (format_data() output)
    \item clump_kb. 10000
    \item clump_r2. 0.001
    \item clump_p1. 1
    \item clump_p2. 1
    \item pop. "EUR" ("EUR", "SAS", "EAS", "AFR", "AMR")
  }
}
\item{harmonise.args}{a list of options appropriate for harmonise_data(),
  \enumerate{
    \item exposure_data. NULL (exposure data from format_data())
    \item outcome_data. NULL (outcome data from format_data())
    \item action. 2 (1, 2, 3)
  }
}
\item{prefix}{a prefix for output files}
\item{...}{other options to be added}
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
# SCALLOP-INF
INF <- Sys.getenv("INF")
outcomes <- "ukb-b-19732"
prot <- "MMP.10"
type <- "cis"
wd <- file.path(INF,"mr",type)
f <- paste0(prot,"-",type,".mrx")
d <- file.path(wd,f)
library(TwoSampleMR)
\dontrun{
run_TwoSampleMR(exposure.args=format_data.args(d),outcome.args=extract_outcome_data.args(outcomes=outcomes),
                clump.args=clump_data.args(),harmonise.args=harmonise_data.args(),
                prefix=file.path(wd,paste(outcomes,prot,type,sep="-")))
}
# pQTLtools
d <- file.path(find.package("pQTLtools"),"tests",f)
run_TwoSampleMR(exposure.args=format_data.args(d),outcome.args=extract_outcome_data.args(outcomes=outcomes),
                clump.args=clump_data.args(),harmonise.args=harmonise_data.args(),
                prefix=paste(outcomes,prot,type,sep="-"))
}
\note{Adapted from script by Dimou NL, Tsilidis KK}
\keyword{utilities}
