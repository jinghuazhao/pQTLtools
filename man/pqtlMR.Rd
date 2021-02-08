\name{pqtlMR}
\alias{pqtlMR}
\title{Bidirectional pQTL-MR analysis}
\usage{pqtlMR(Ins,Ids,prefix="INF1",reverse=FALSE)}
\description{
Given a list of pQTLs as instrumental variables (IVs), this function
performs MR analysis according to epigraphdb.
}
\arguments{
\item{Ins}{a data.frame containing a list of pQTLs and their effect sizes}
\item{Ids}{a vector containing a list of MR-Base IDs}
\item{prefix}{a prefix for output files}
\item{reverse}{if TRUE, perform reverse MR}
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
ivs <- read.table(ms,as.is=TRUE,header=TRUE)
ids <- "ieu-b-18"
# MR
pqtlMR(ivs,ids,prefix="MS")
# reverse MR
pqtlMR(ivs,ids,prefix="MS_rev",reverse=TRUE)
dir()
# prefix=INF1 for the SCALLOP/INF1
prefix <- Sys.getenv("prefix")
ivs <- within(read.table(paste0(prefix,".ins"),as.is=TRUE,header=TRUE),{pval=10^pval})
ids <- scan("ieu-id.txt",what="")
pQTLtools::pqtlMR(ivs,ids)
}
}
\note{adapted from script by Jie Zheng}
\keyword{utilities}
