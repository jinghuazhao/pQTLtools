\name{make_ExpressionSet}
\alias{make_ExpressionSet}
\title{A call to expressionSet class}
\usage{make_ExpressionSet(exprs,phenoData,experimentData=NULL,annotation=NULL)}
\description{
This is really a direct call to the Bioconductor/Biobase class.
}
\arguments{
\item{exprs}{Expression data}
\item{phenoData}{Phenotype}
\item{experimentData}{Information on data source}
\item{annotation}{Annotation information}
}
\details{
The explicit call make it easier to handle proteomic data for other downstream analyses.
}

\value{
An ExpressionSet object.
}

\examples{
dataDirectory <- system.file("extdata", package="Biobase")
exprsFile <- file.path(dataDirectory, "exprsData.txt")
exprs <- as.matrix(read.table(exprsFile, header=TRUE, sep="\t", row.names=1, as.is=TRUE))
pDataFile <- file.path(dataDirectory, "pData.txt")
pData <- read.table(pDataFile, row.names=1, header=TRUE, sep="\t")
all(rownames(pData)==colnames(exprs))
metadata <- data.frame(labelDescription=
                       c("Patient gender",
                         "Case/control status",
                         "Tumor progress on XYZ scale"),
                       row.names=c("gender", "type", "score"))
phenoData <- new("AnnotatedDataFrame", data=pData, varMetadata=metadata)
experimentData <- new("MIAME",
  name="Pierre Fermat",
  lab="Francis Galton Lab",
  contact="pfermat@lab.not.exist",
  title="Smoking-Cancer Experiment",
  abstract="An example ExpressionSet",
  url="www.lab.not.exist",
  other=list(notes="Created from text files"))
exampleSet <- make_ExpressionSet(exprs,phenoData,experimentData,annotation="hgu95av2")
require(Biobase)
data(sample.ExpressionSet)
identical(exampleSet,sample.ExpressionSet)
esApply(exampleSet,1,hist)
lm(score~gender+X31739_at,data=exampleSet)
}
\note{Adapted from Bioconductor/Biobase following a number of proteomic pilot studies}
\keyword{utilities}
