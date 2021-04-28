\name{pqtlMR}
\alias{pqtlMR}
\title{Bidirectional pQTL-MR analysis}
\usage{pqtlMR(pqtlMRinput, plot=FALSE, prefix="pQTL-combined-",reverse=FALSE)}
\description{
Given a list of pQTLs as instrumental variables (IVs), this function performs MR analysis according to epigraphdb.
}
\arguments{
\item{pqtlMRinput}{A list of the following elements for greater flexibility,
  \enumerate{
    \item Ins result of format_data()
    \item Ids result of extract_outcome_data()
    \item harmonise result of harmonise_data()
  }
}
\item{plot}{to produce plots}
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
Zheng J, et al. (2020). Phenome-wide Mendelian randomization mapping the influence of the plasma proteome on complex diseases. *Nature Genetics* 52(10): 1122-1131.
}

\examples{
library(TwoSampleMR)
library(pQTLtools)
# Original examples
ins <- read.csv(file.path(find.package("pQTLtools",lib.loc=.libPaths()),"tests","Ins.csv"))
ids <- c("ieu-a-7","ebi-a-GCST007432")
Ins <- format_data(ins)
Ids <- extract_outcome_data(snps=with(Ins,SNP),outcomes=ids)
harmonise <- harmonise_data(Ins,Ids)
input <- list(Ins, Ids, harmonise)
pqtlMR(input,plot=FALSE)
# A bidirectional analysis
ms <- read.table(file.path(find.package("pQTLtools",lib.loc=.libPaths()),"tests","ms.ins"),
                 header=TRUE)
Ins <- format_data(ms,samplesize_col="N")
Ids <- extract_outcome_data(snps=with(Ins,SNP),outcomes="ieu-b-18")
harmonise <- harmonise_data(Ins,Ids)
input <- list(Ins, Ids, harmonise)
# MR
pqtlMR(input,prefix="MS-")
# reverse MR
pqtlMR(input,prefix="MS_rev-",reverse=TRUE)
}
\note{Adapted from script by Jie Zheng}
\keyword{utilities}
