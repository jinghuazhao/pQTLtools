# A call to expressionSet class

This is really a direct call to the Bioconductor/Biobase class.

## Usage

``` r
make_ExpressionSet(
  assayData,
  phenoData = Biobase::annotatedDataFrameFrom(assayData, byrow = FALSE),
  featureData = Biobase::annotatedDataFrameFrom(assayData, byrow = TRUE),
  experimentData = Biobase::MIAME(),
  annotation = character(),
  protocolData = Biobase::annotatedDataFrameFrom(assayData, byrow = FALSE),
  ...
)
```

## Arguments

-   assayData:

    Expression data.

-   phenoData:

    Phenotype.

-   featureData:

    featureData.

-   experimentData:

    Information on data source.

-   annotation:

    Annotation information.

-   protocolData:

    protocol information.

-   ...:

    Other options.

## Value

An ExpressionSet object.

## Details

The explicit call make it easier to handle proteomic data for other
downstream analyses.

## Note

Adapted from Bioconductor/Biobase following a number of proteomic pilot
studies.

## Examples

``` r
dataDirectory <- system.file("extdata", package="Biobase")
exprsFile <- file.path(dataDirectory, "exprsData.txt")
exprs <- as.matrix(read.table(exprsFile, header=TRUE, sep="\t", row.names=1, as.is=TRUE))
pDataFile <- file.path(dataDirectory, "pData.txt")
pData <- read.table(pDataFile, row.names=1, header=TRUE, sep="\t")
all(rownames(pData)==colnames(exprs))
#> [1] TRUE
metadata <- data.frame(labelDescription=
                       c("Patient gender",
                         "Case/control status",
                         "Tumor progress on XYZ scale"),
                       row.names=c("gender", "type", "score"))
suppressMessages(library(Biobase))
phenoData <- Biobase::AnnotatedDataFrame(data=pData, varMetadata=metadata)
experimentData <- Biobase::MIAME(
    name="Pierre Fermat",
    lab="Francis Galton Lab",
    contact="pfermat@lab.not.exist",
    title="Smoking-Cancer Experiment",
    abstract="An example ExpressionSet",
    url="www.lab.not.exist",
    other=list(notes="Created from text files"))
exampleSet <- pQTLtools::make_ExpressionSet(exprs,phenoData,
                                 experimentData=experimentData,
                                 annotation="hgu95av2")
data(sample.ExpressionSet, package="Biobase")
identical(exampleSet,sample.ExpressionSet)
#> [1] FALSE
invisible(Biobase::esApply(exampleSet,2,hist))


























lm(score~gender+X31739_at,data=exampleSet)
#> 
#> Call:
#> lm(formula = score ~ gender + X31739_at, data = exampleSet)
#> 
#> Coefficients:
#> (Intercept)   genderMale    X31739_at  
#>   0.6006673    0.0108515   -0.0003012  
#> 
```
