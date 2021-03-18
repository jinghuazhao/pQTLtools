\name{pqtlMR}
\alias{pqtlMR}
\title{Bidirectional pQTL-MR analysis}
\usage{pqtlMR(Ins=format_data.args(),Ids=extract_outcome_data.args(),prefix="INF1",reverse=FALSE)}
\description{
Given a list of pQTLs as instrumental variables (IVs), this function performs MR analysis according to epigraphdb.
}
\arguments{
\item{Ins}{a list of options appropriate for format_data(),
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
\item{Ids}{a list of options appropriate for extract_outcome_data(),
  \enumerate{
    \item snp. snpname
    \item outcomes. outcome id
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
\item{prefix}{a prefix for output files}
\item{reverse}{if TRUE, perform reverse MR}
\item{...}{other possible options}
}
\details{
This function is based on TwoSampleMR.
}

\value{
No value is returned but several files.
}

\section{References}{
Zheng, J., et al. (2020). Phenome-wide Mendelian randomization mapping the influence of the plasma proteome on complex diseases. *Nature Genetics* 52(10): 1122-1131.
}

\examples{
\dontrun{
# A bidirectional MR analysis
ms <- file.path(find.package("pQTLtools"),"tests","ms.ins")
ivs <- format_data.args(file=ms)
ids <- extract_outcome_data.args(outcomes="ieu-b-18")
# MR
pqtlMR(ivs,ids,prefix="MS")
# reverse MR
pqtlMR(ivs,ids,prefix="MS_rev",reverse=TRUE)
dir()
# prefix=INF1 for the SCALLOP/INF1
prefix <- Sys.getenv("prefix")
ivs <- format_data.args(file=paste0(prefix,".ins",log_pval=TRUE))
ids <- extract_outcome_data.args(outcomes=scan("ieu-id.txt",what=""))
pQTLtools::pqtlMR(ivs,ids)
}
}
\note{Adapted from script by Jie Zheng}
\keyword{utilities}
