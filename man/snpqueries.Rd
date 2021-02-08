\name{snpqueries}
\alias{snpqueries}
\title{phenoscanner snpqueries in batches}
\usage{snpqueries(snplist,catalogue="pQTL",proxies="EUR",p=5e-8,r2=0.8,build=37,wait=TRUE)}
\description{
R/phenoscanner only allows for certain number of items supplied. This simple function return
a large number of calls in batches as well as generating SNPIDs.
}
\arguments{
\item{snplist}{a list of SNPs}
\item{catalogue}{"None","eQTL","mQTL","methQTL","pQTL","GWAS"}
\item{proxies}{"None", "AFR","AMR","EAS","EUR","SAS"}
\item{p}{p value threshold}
\item{r2}{r2 for LD}
\item{build}{37, 38}
\item{wait}{a flag to wait for 1hr for every 500 SNPs}
}
\details{
Batches are generated and queries are combined into one.
}

\value{
The returned value is a list containing snps and results:
}

\section{References}{
Sun BB, et al. (2018). Genomic atlas of the human plasma proteome. *Nature* 558: 73-79.
}
\seealso{
\code{\link[phenoscanner]{phenoscanner}}
}

\examples{
\dontrun{
# single SNP
snpqueries("rs704")

# SCALLOP/INF
INF <- Sys.getenv("INF")
rsid <- scan(paste(INF,'work','INF1.merge.snp',sep='/'),"")
r <- snpqueries(rsid,catalogue='pQTL',p=1e-11)
INTERVAL_Olink <- subset(with(r,results),efo=='EFO_0004747' & pmid=='29875488')
save(INTERVAL_Olink,file='INTERVAL_Olink.rda',compress='xz')

# --- query intersect proteins ---

# SomaLogic intersect
SomaLogic_overlap_list <- subset(st4,UniProt \%in\% intersect_list)
r <- snpqueries(SomaLogic_overlap_list[,6],catalogue='pQTL',p=1e-11)
SomaLogic_overlap <- subset(with(r,results),efo=='EFO_0004747' & pmid=='29875488')
save(SomaLogic_overlap_list,SomaLogic_overlap,file='SomaLogic_overlap.rda',compress='xz')
SomaLogic_result <- merge(SomaLogic_overlap_list,SomaLogic_overlap,
                          by.x="Sentinel.variant*",by.y="snp")

# Olink intersect
INF1_merge_rsid <- read.delim(paste(INF,"work","INF1.merge-rsid",sep="/"))
INF1_merge_rsid_uniprot <- merge(INF1_merge_rsid,inf1,by="prot")
Olink_overlap_list <- subset(INF1_merge_rsid_uniprot,uniprot \%in\% intersect_list)
r <- snpqueries(with(Olink_overlap_list,MarkerName),catalogue='pQTL',p=1e-11)
Olink_overlap <- subset(with(r,results),efo=='EFO_0004747' & pmid=='29875488')
save(Olink_overlap_list,Olink_overlap,file='Olink_overlap.rda',compress='xz')
Olink_result <- merge(Olink_overlap_list,Olink_overlap,by.x="MarkerName",by.y="snp")
}
}
\note{adapted from custom codings}
\keyword{utilities}
