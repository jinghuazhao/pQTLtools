\name{uniprot2ids}
\alias{uniprot2ids}
\title{UniProt IDs to others}
\usage{uniprot2ids(uniprotid,to,query)}
\description{
A function which converts UniProt IDs to others.
}
\arguments{
\item{uniprotid}{Source IDs}
\item{to}{To IDs}
\item{query}{A query}
}
\details{
This function is based on the Python3 script from UniProt.
}

\value{
A UniProt-ID mapping
}

\section{References}{
See https://www.uniprot.org/help/api_idmapping
}

\examples{
\dontrun{
uniprotid <- "ACC+ID"
to <- "CHEMBL_ID"
query <- noquote(inf1[["uniprot"]])
query <- paste(query,collapse=" ")
r <- pQTLtools::uniprot2ids(uniprotid,to,query)
cat(r,file="INF1.merge.chembl")
}
}
\note{Adapted from script by UniProt}
\keyword{utilities}
