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
#' \dontrun{
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
#' phenoData <- Biobase::AnnotatedDataFrame(data=pData, varMetadata=metadata)
#' experimentData <- Biobase::MIAME(
#'     name="Pierre Fermat",
#'     lab="Francis Galton Lab",
#'     contact="pfermat@lab.not.exist",
#'     title="Smoking-Cancer Experiment",
#'     abstract="An example ExpressionSet",
#'     url="www.lab.not.exist",
#'     other=list(notes="Created from text files"))
#' exampleSet <- pQTLtools::make_ExpressionSet(exprs,phenoData,
#'                                  experimentData=experimentData,
#'                                  annotation="hgu95av2")
#' data(sample.ExpressionSet, package="Biobase")
#' identical(exampleSet,sample.ExpressionSet)
#' invisible(Biobase::esApply(exampleSet,2,hist))
#' lm(score~gender+X31739_at,data=exampleSet)
#' }
#' @note
#' Adapted from Bioconductor/Biobase following a number of proteomic pilot studies.
#' @keywords utilities

make_ExpressionSet <- function(assayData,
                      phenoData=Biobase::annotatedDataFrameFrom(assayData, byrow=FALSE),
                      featureData=Biobase::annotatedDataFrameFrom(assayData, byrow=TRUE),
                      experimentData=Biobase::MIAME(),
                      annotation=character(),
                      protocolData=Biobase::annotatedDataFrameFrom(assayData, byrow=FALSE),...)
{
  Biobase::ExpressionSet(assayData,phenoData=phenoData,
                         featureData=featureData,
                         experimentData=experimentData,
                         annotation=annotation,
                         protocolData=protocolData,...)
}

#' Limit of detection analysis for ExpressionSet objects
#'
#' Computes the percentage of values below a feature-specific lower limit of detection (LOD)
#' in an \code{ExpressionSet}. Commonly used in proteomic quality assessment.
#'
#' @param eset An \code{ExpressionSet} object containing expression values and feature metadata.
#' Must include a numeric column \code{lod.max} in \code{fData(eset)}.
#'
#' @param flagged Character string indicating whether flagged samples should be removed.
#' One of \code{"OUT"} (remove flagged samples) or \code{"IN"} (retain all samples).
#'
#' @return An \code{ExpressionSet} with an added feature annotation column:
#' \item{pc.belowLOD.new}{Percentage of values below the LOD per feature.}
#'
#' @details
#' The function compares expression values from \code{exprs(eset)} to feature-specific
#' LOD thresholds stored in \code{fData(eset)$lod.max}. Computation is vectorised using
#' \code{sweep()} for efficiency and consistency with Bioconductor conventions.
#'
#' Flagged samples (if present in \code{pData(eset)$Flagged}) can optionally be removed.
#'
#' @examples
#' \dontrun{
#' suppressMessages(library(Biobase))
#' data(sample.ExpressionSet, package = "Biobase")
#' exampleSet <- sample.ExpressionSet
#' set.seed(1)
#' probs <- runif(nrow(exampleSet))
#' Biobase::fData(exampleSet)$lod.max <- vapply(seq_len(nrow(exampleSet)), function(i)
#'            quantile(Biobase::exprs(exampleSet)[i, ], probs = probs[i], na.rm = TRUE),
#'            numeric(1))
#' lod <- get.prop.below.LLOD(exampleSet)
#' fd <- Biobase::fData(lod)
#' fd <- fd[order(fd$pc.belowLOD.new, decreasing = TRUE), ]
#' knitr::kable(head(fd))
#' plot(fd$pc.belowLOD.new, main = "Random quantile cut off", ylab = "% below LOD")
#' }
#' @export
#'
get.prop.below.LLOD <- function(eset, flagged = c("OUT", "IN")) {
    flagged <- match.arg(flagged)
    if (!inherits(eset, "ExpressionSet")) {
        stop("'eset' must be an ExpressionSet object")
    }
    fd <- Biobase::fData(eset)
    pd <- Biobase::pData(eset)
    expr_mat <- Biobase::exprs(eset)
    # check required column
    if (!"lod.max" %in% colnames(fd)) {
        stop("Feature data must contain column 'lod.max'")
    }
    lod <- fd$lod.max
    # optional removal of flagged samples
    if ("Flagged" %in% colnames(pd)) {
        ind_fl <- which(pd$Flagged == "Flagged")

        if (flagged == "OUT" && length(ind_fl) > 0) {
            eset <- eset[, -ind_fl]
            expr_mat <- Biobase::exprs(eset)
        }
    }
    fd <- Biobase::fData(eset)
    lod <- fd$lod.max
    # dimension safety check
    if (nrow(expr_mat) != length(lod)) {
        stop("Mismatch: number of features in exprs does not match length of lod.max")
    }
    # vectorised core computation (Bioconductor standard)
    below <- sweep(expr_mat, 1, lod, FUN = "<=")
    denom <- rowSums(!is.na(expr_mat))
    denom[denom == 0] <- NA_real_
    fd$pc.belowLOD.new <- 100 * rowSums(below, na.rm = TRUE) / denom
    Biobase::fData(eset) <- fd
    eset
}
