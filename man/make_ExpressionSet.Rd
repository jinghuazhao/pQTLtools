\name{make_ExpressionSet}
\alias{make_ExpressionSet}
\title{A call to expressionSet class}
\usage{make_ExpressionSet(assayData,
            phenoData=Biobase::annotatedDataFrameFrom(assayData, byrow=FALSE),
            featureData=Biobase::annotatedDataFrameFrom(assayData, byrow=TRUE),
            experimentData=Biobase::MIAME(),
            annotation=character(),
            protocolData=Biobase::annotatedDataFrameFrom(assayData, byrow=FALSE),
           ...
           )
      }
\description{
This is really a direct call to the Bioconductor/Biobase class.
}
\arguments{
\item{assayData}{Expression data}
\item{phenoData}{Phenotype}
\item{featureData}{featureData}
\item{experimentData}{Information on data source}
\item{annotation}{Annotation information}
\item{protocolData}{protocol information}
\item{...}{Other options}
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
require(Biobase)
require(pQTLtools)
phenoData <- new("AnnotatedDataFrame", data=pData, varMetadata=metadata)
experimentData <- new("MIAME",
  name="Pierre Fermat",
  lab="Francis Galton Lab",
  contact="pfermat@lab.not.exist",
  title="Smoking-Cancer Experiment",
  abstract="An example ExpressionSet",
  url="www.lab.not.exist",
  other=list(notes="Created from text files"))
exampleSet <- make_ExpressionSet(exprs,phenoData,experimentData=experimentData,
                                 annotation="hgu95av2")
data(sample.ExpressionSet)
identical(exampleSet,sample.ExpressionSet)
invisible(esApply(exampleSet,2,hist))
lm(score~gender+X31739_at,data=exampleSet)
}
\note{Adapted from Bioconductor/Biobase following a number of proteomic pilot studies}
\keyword{utilities}
