\name{import_OpenGWAS}
\alias{import_OpenGWAS}
\title{Import OpenGWAS}
\usage{import_OpenGWAS(opengwas_id, region, verbose = TRUE)}
\description{
A function which imports OpenGWAS.
}
\arguments{
\item{opengwas_id}{An OpenGWAS id}
\item{region}{chr:start-end}
\item{verbose}{Extra information}
}
\details{
This function is derived from SCALLOP/INF work.
}

\value{
A summary statistic object
}

\section{References}{
Lyon M, Andrews SJ, Elsworth B, Gaunt TR, Hemani G, Marcora E. The variant call format provides efficient and robust storage of GWAS summary statistics. bioRxiv 2020.05.29.115824; doi: https://doi.org/10.1101/2020.05.29.115824
}

\examples{
\dontrun{
options(width=200)
gwasvcf::set_bcftools(path=file.path(HPC_WORK,"bin","bcftools"))
# MPV ARHGEF3 region
opengwas_id <- "ebi-a-GCST004599"
region <- "3:56649749-57049749"
mpv_ARHGEF3 <- import_OpenGWAS(opengwas_id,region)
# all immune-related
INF <- Sys.getenv("INF")
HPC_WORK <- Sys.getenv("HPC_WORK")
opengwas_ids <- scan(file.path(INF,"OpenGWAS","ieu.list"),what="")
unavail <-c("ieu-b-18","finn-a-M13_SLE","finn-a-D3_SARCOIDOSIS")
opengwas_ids <- subset(opengwas_ids,!opengwas_ids\%in\%unavail)
region <- "1:100-2000000"
library(pQTLtools)
summary_list = purrr::map(opengwas_ids[1:2], ~import_OpenGWAS(., region))
}
}
\note{Adapted function}
\keyword{utilities}
