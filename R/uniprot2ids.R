#' UniProt IDs to others
#'
#' A function which converts UniProt IDs to others.
#'
#' @param uniprotid Source IDs.
#' @param to To IDs.
#' @param query A query.
#'
#' @details
#' This function is based on the Python3 script from UniProt.
#' See https://www.uniprot.org/help/api_idmapping
#'
#' @return A UniProt-ID mapping
#'
#' @examples
#' \dontrun{
#' uniprotid <- "ACC+ID"
#' to <- "CHEMBL_ID"
#' query <- noquote(gap.datasets::inf1[["uniprot"]])
#' query <- paste(query,collapse=" ")
#' r <- pQTLtools::uniprot2ids(uniprotid,to,query)
#' cat(r,file="INF1.merge.chembl")
#' }
#' @note
#' Adapted from script by UniProt
#' @keywords utilities

uniprot2ids <- function(uniprotid="ACC+ID",to,query)
{
  rt <- find.package("pQTLtools")
  f <- file.path(rt ,"UniProt","uniprot2ids.py")
  reticulate::source_python(f)
  invisible(uniprot2ids(uniprotid,to,query))
}
