#' A call to expressionSet class
#'
#' This is really a direct call to the Bioconductor/Biobase class.
#'
#' @param assayData Expression data.
#' @param phenoData Phenotype.
#' @param featureData featureData.
#' @param experimentData Information on data source.
#' @param annotation Annotation information.
#' @param protocolData protocol information.
#' @param ... Other options.
#'
#' @details
#' The explicit call make it easier to handle proteomic data for other downstream analyses.
#'
#' @export
#' @return
#' An ExpressionSet object.
#'
#' @examples
#' dataDirectory <- system.file("extdata", package="Biobase")
#' exprsFile <- file.path(dataDirectory, "exprsData.txt")
#' exprs <- as.matrix(read.table(exprsFile, header=TRUE, sep="\t", row.names=1, as.is=TRUE))
#' pDataFile <- file.path(dataDirectory, "pData.txt")
#' pData <- read.table(pDataFile, row.names=1, header=TRUE, sep="\t")
#' all(rownames(pData)==colnames(exprs))
#' metadata <- data.frame(labelDescription=
#'                        c("Patient gender",
#'                          "Case/control status",
#'                          "Tumor progress on XYZ scale"),
#'                        row.names=c("gender", "type", "score"))
#' suppressMessages(library(Biobase))
#' suppressMessages(library(pQTLtools))
#' phenoData <- new("AnnotatedDataFrame", data=pData, varMetadata=metadata)
#' experimentData <- new("MIAME",
#'   name="Pierre Fermat",
#'   lab="Francis Galton Lab",
#'   contact="pfermat@lab.not.exist",
#'   title="Smoking-Cancer Experiment",
#'   abstract="An example ExpressionSet",
#'   url="www.lab.not.exist",
#'   other=list(notes="Created from text files"))
#' exampleSet <- make_ExpressionSet(exprs,phenoData,experimentData=experimentData,
#'                                  annotation="hgu95av2")
#' data(sample.ExpressionSet)
#' identical(exampleSet,sample.ExpressionSet)
#' invisible(esApply(exampleSet,2,hist))
#' lm(score~gender+X31739_at,data=exampleSet)
#'
#' @note
#' Adapted from Bioconductor/Biobase following a number of proteomic pilot studies.
#' @keywords utilities

make_ExpressionSet <- function(assayData,
                      phenoData=Biobase::annotatedDataFrameFrom(assayData, byrow=FALSE),
                      featureData=Biobase::annotatedDataFrameFrom(assayData, byrow=TRUE),
                      experimentData=Biobase::MIAME(),
                      annotation=character(),
                      protocolData=Biobase::annotatedDataFrameFrom(assayData, byrow=FALSE),...)
Biobase::ExpressionSet(assayData,phenoData=phenoData,
                       featureData=featureData,
                       experimentData=experimentData,
                       annotation=annotation,
                       protocolData=protocolData,...)

#' Limit of detection analysis
#'
#' The function obtains lower limit of detection as in proteomic analysis.
#'
#' @param eset An ExpressionSet object.
#' @param flagged A flag is an indicator for sample exclusion.
#' @export
#' @return An updated ExpressionSet object.
#' @examples
#' suppressMessages(library(Biobase))
#' data(sample.ExpressionSet)
#' exampleSet <- sample.ExpressionSet
#' fData(exampleSet)
#' fData(exampleSet)$lod.max <- apply(exprs(exampleSet),1,quantile,runif(nrow(exampleSet)))
#' lod <- get.prop.below.LLOD(exampleSet)
#' x <- dplyr::arrange(fData(lod),desc(pc.belowLOD.new))
#' knitr::kable(head(lod))
#' plot(x[,2], main="Random quantile cut off", ylab="<lod%")
#' @author James Peters

get.prop.below.LLOD <- function(eset, flagged = 'OUT'){

  ## A function to calculate no. of proteins i.e NA per sample (missing or <LLD per sample)
  # arguments 'eset' and 'flagged'
  # flagged indicates whether Flagged samples should be excluded (if they have not been already)

  if(!inherits(eset, what= "ExpressionSet")){
    stop("'eset' argument must inherit class ExpressionSet")
  }

  if (!flagged %in% c('IN','OUT')){
    stop("'flagged' argument must be 'IN' or 'OUT")
  }

  for (p in "stringr") {
    if (length(grep(paste("^package:", p, "$", sep=""), search())) == 0) {
       if (!requireNamespace(p, quietly = TRUE))
           warning(paste("This function needs package `", p, "' to be fully functional; please install", sep=""))
    }
  }

  # best to cut flagged samples first at eset stage:
  # risk of messing up if cutting from matrix, and then dont edit pData

  ind.fl <- which(eset$Flagged == 'Flagged')

  if (flagged == "IN"){

    if (length(ind.fl) > 0){
      panel <- unique(eset$panel)
      mytit <- paste(toupper(panel), "panel \n (flagged retained)")
    } else{
      panel <- unique(eset$panel)
      mytit <- paste(toupper(panel), "panel \n (no. flagged samples = 0)")
    }

  } else if (flagged == "OUT"){
    # cut flagged samples

    if (length(ind.fl) > 0){
      eset <- eset[, -ind.fl] # nb annoying ESet behaviour: cols are samples
      panel <- unique(eset$panel)
      mytit <- paste(toupper(panel), "panel \n (flagged removed)")
    } else{
      panel <- unique(eset$panel)
      mytit <- paste(toupper(panel), "panel \n (no. flagged samples = 0)")
    }

  }

  E <- t(Biobase::exprs(eset))

  p.annot <- Biobase::fData(eset)

  p.annot$pc.belowLOD.new <- NA

  # % proteins in each sample
  ##miss.by.prot <- apply(E, 2, FUN=function(x) 100*sum(is.na(x))/nrow(E) )

  for (i in 1:ncol(E)){
    m <- sum( E[,i] <= p.annot$lod.max[i], na.rm=T ) # number of samples <= LLOD
    t <- length( na.omit(E[,i] )) # denominator NB use this rather than just nrow(E) to make code robust in event of missing values
    p.annot$pc.belowLOD.new[i] <- 100*m/t
  }

  Biobase::fData(eset) <- p.annot

  eset
  #eof
}
