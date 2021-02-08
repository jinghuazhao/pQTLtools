\name{pQTLtools-internal}
\title{Internal data/functions for pQTLtools}
\alias{biomaRt}
\alias{caprion}
\alias{hg19}
\alias{hg19Tables}
\alias{inf1}
\alias{Olink_NGS}
\alias{Olink_qPCR}
\alias{run_coloc}
\alias{SomaLogic160410}
\alias{st4}
\alias{st6}
\alias{swath_ms}
\description{
A list of protein panels: caprion, hg10Tables, Olink_NGS, Olink_qPCR (inf1 is heavily edited), SomaLogic160410, swath_ms.

In cases of incomplete information, datasets called hg19Tables based on UCSC and biomaRt from Ensembl are provided.
}
\arguments{
  \item{chr}{chromosome}
  \item{start}{Start position}
  \item{end}{End position}
  \item{gene}{Gene name}
  \item{UniProt}{UniProt ID}
}
\examples{
\dontrun{
# datasets
head(biomaRt)

# Olink-SomaLogic panel overlap
p <- list(setdiff(inf1$uniprot,"P23560"),
          setdiff(SomaLogic160410$UniProt[!is.na(SomaLogic160410$UniProt)],"P23560"))
cnames <- c("INF1","SomaLogic")
VennDiagram::venn.diagram(x = p, category.names=cnames,
                          filename='os.png', imagetype="png", output=TRUE)
m <- merge(inf1,SomaLogic160410,by.x="uniprot",by.y="UniProt")
u <- setdiff(with(m,unique(uniprot)),"P23560")
options(width=220)
o <- subset(inf1,uniprot \%in\% u)
dim(o)
o
vars <- c("UniProt","chr","start","end","extGene","Target","TargetFullName")
s <- subset(SomaLogic160410[vars], UniProt \%in\% u)
dim(s)
us <- s[!duplicated(s),]
dim(us)
us

# SCALLOP/INF1
INF <- Sys.getenv("INF")
INF1_merge <- merge(inf1,
                    read.delim(file.path(INF,"work","INF1.merge-rsid"),as.is=TRUE),
                    by="prot")
INF1_uniprot <- unique(with(INF1_merge,uniprot))

# INTERVAL SomaLogic at box
HOME <- Sys.getenv("HOME")
box <- read.delim(file.path(HOME,"SomaLogic","doc","INTERVAL-box.tsv"),as.is=TRUE)
box_INF1 <- subset(box,UniProt \%in\% INF1_uniprot)
box_uniprot <- setdiff(unique(with(box_INF1,UniProt)),"P23560")
setdiff(INF1_uniprot,box_uniprot)

# Phenoscanner database
ps <- merge(subset(read.delim(file.path(INF,"work","pQTL_2018.txt.gz"),as.is=TRUE),pmid==29875488),
            box,by.x="trait",by.y="TargetFullName")
z <- subset(ps,UniProts \%in\% INF1_uniprot & p<=1.5e-11)

# ST4 at Nature
st4regions <- subset(st4, UniProt \%in\% INF1_uniprot)
unique_uniprot_list <- setdiff(intersect(st4$UniProt,inf1$uniprot),"P23560")
subset(INF1_merge,uniprot \%in\% unique_uniprot_list)
}
}
\details{This seeds for prototype for extensions}
\keyword{internal}
