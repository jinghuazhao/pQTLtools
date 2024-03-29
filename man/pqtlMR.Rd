% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{pqtlMR}
\alias{pqtlMR}
\title{Basic pQTL-MR analysis}
\usage{
pqtlMR(ivs, ids, mr_plot = FALSE, prefix = "pQTL-combined-", reverse = FALSE)
}
\arguments{
\item{ivs}{Instrumental variables from format_data().}

\item{ids}{A list of MR-Base IDs.}

\item{mr_plot}{to produce plots.}

\item{prefix}{a prefix for output files.}

\item{reverse}{if TRUE, perform reverse MR.}
}
\value{
No value is returned but several files.
}
\description{
This function takes data intrumental variables as produced by \code{format_data()} and
used to perform MR analysis against a list of outcomes from MR-Base.
}
\details{
This function is adapted from \insertCite{zheng20;textual}{pQTLtools} and in
spirit similar to \code{run_TwoSampleMR}.
}
\note{
Adapted from script by Jie Zheng.
}
\examples{
library(TwoSampleMR)
library(pQTLtools)
# Original examples
f <- file.path(find.package("pQTLtools",lib.loc=.libPaths()),"tests","Ins.csv")
ivs <- format_data(read.csv(f))
ids <- c("ieu-a-7","ebi-a-GCST007432")
pqtlMR(ivs, ids, mr_plot=FALSE)
# A bidirectional analysis
f <- file.path(find.package("pQTLtools",lib.loc=.libPaths()),"tests","ms.ins")
ivs <- format_data(read.table(f, header=TRUE), samplesize_col="N")
ids <- "ieu-b-18"
# MR
pqtlMR(ivs, ids, prefix="MS-")
# reverse MR
pqtlMR(ivs, ids, ,prefix="MS_rev-",reverse=TRUE)

}
\references{
\insertAllCited{}
}
\keyword{utilities}
