#' Variant consequence
#'
#' This function maps consequences (CSQ) of a given list of variants to those in an annotated
#' set based on facilities such as variant effect predictor (VEP). In the case of
#' protein-altering variants (PAVs), many consequences can be involved. The procedure also
#' takes into account linkage disequilibrium (LD) in flanking regions.
#'
#' @param query_loci A data.frame of loci whose consequences are to be obtained.
#' @param annotated_loci A data.frame of annotated loci.
#' @param pattern A character string or pattern to match the consequences against.
#' @param ldops Arguments for `ieugwasr::ld_matrix_local()`, typically a list with `bfile` and `plink` paths.
#' @param flanking A numeric value specifying the flanking distance (default is 1e6).
#' @param pop A character string specifying the reference population for `ieugwasr::ld_matrix()` (default is "EUR").
#' @param verbose A logical flag to show nonexistent variants (default is TRUE).
#'
#' @return A data.frame with an indicator showing whether the consequences are present or absent.
#' @export
#'
#' @examples
#' \dontrun{
#' options(width=2000)
#' suppressMessages(require(dplyr))
#' suppressMessages(require(stringr))
#' # SCALLOP-INF list
#' METAL <- read.delim(file.path(find.package("pQTLtools"), "tests", "INF1.METAL")) %>%
#'          dplyr::left_join(gap.datasets::inf1[c("prot", "gene")]) %>%
#'          dplyr::mutate(prot = gene, chr = Chromosome, pos = Position)
#' # VEP output
#' vep <- "/rds/project/rds-zuZwCZMsS0w/Caprion_proteomics/analysis/bgen/vep"
#' pattern <- paste(protein_altering_variants, collapse = "|")
#' suppressMessages(require(GenomicRanges))
#' INF <- "/rds/project/rds-zuZwCZMsS0w/olink_proteomics/scallop/INF"
#' plink <- "/rds/user/jhz22/hpc-work/bin/plink"
#' # CSQ
#' b <- list()
#' for (i in unique(dplyr::pull(METAL, Chromosome))) {
#'   m <- dplyr::filter(METAL, Chromosome %in% i) %>%
#'        dplyr::select(chr, pos, MarkerName, prot) %>%
#'        dplyr::mutate(rsid = gsub("chr", "", MarkerName)) %>%
#'        dplyr::select(-MarkerName)
#'   u <- read.delim(file.path(vep, paste0("chr", i, ".tab.gz"))) %>%
#'        dplyr::select(Chrom, Pos, X.Uploaded_variation, Consequence) %>%
#'        setNames(c("chr", "pos", "rsid", "csq"))
#'   bfile <- file.path(INF, "INTERVAL", "per_chr", paste0("snpid", i))
#'   b[[i]] <- csq(m, u, pattern, ldops = list(bfile = bfile, plink = plink))
#' }
#' r <- dplyr::bind_rows(b) %>%
#'      dplyr::filter(r2 >= 0.8) %>%
#'      dplyr::rename(gene = prot) %>%
#'      dplyr::mutate(seqnames = as.integer(seqnames), pos = as.integer(pos)) %>%
#'      dplyr::arrange(seqnames, pos) %>%
#'      dplyr::select(-ref.seqnames, -ref.start, -ref.end, -seqnames, -pos)
#' }
#'
csq <- function(query_loci, annotated_loci, pattern, ldops=NULL, flanking=1e6, pop="EUR", verbose=TRUE) {
  rsid <- seqnames <- start <- strand <- width <- bfile <- plink <- NA
  
  subject <- with(query_loci,
    GenomicRanges::GRanges(seqnames=chr, IRanges::IRanges(start=pos-flanking, end=pos+flanking),
                           rsid=rsid, pos=pos, prot=prot))
  
  query <- with(filter(annotated_loci, stringr::str_detect(csq, pattern)),
    GenomicRanges::GRanges(seqnames=chr, IRanges::IRanges(start=pos-flanking, end=pos+flanking),
                           rsid=rsid, pos=pos, csq=csq))
  
  ov <- GenomicRanges::findOverlaps(query, subject) %>%
    data.frame()
  
  ov1 <- subject[ov$subjectHits, ] %>%
    data.frame() %>%
    dplyr::select(-strand, -width)
  ov2 <- query[ov$queryHits, ] %>%
    data.frame() %>%
    dplyr::select(-strand, -width) %>%
    dplyr::mutate(rsid=if_else(rsid=="-", paste0("chr", seqnames, ":", start), rsid))
  
  b <- dplyr::bind_cols(data.frame(ov1),
                        data.frame(ov2) %>% setNames(paste("ref", names(ov2), sep=".")))
  b_granges <- GenomicRanges::makeGRangesFromDataFrame(b, keep.extra.columns = TRUE)
  variant_list <- unique(c(b[["rsid"]], b[["ref.rsid"]]))
  
  if (!is.null(ldops)) {
    r <- ieugwasr::ld_matrix_local(variants = c(b[["rsid"]], b[["ref.rsid"]]),
                                   bfile = ldops[["bfile"]],
                                   plink_bin = ldops[["plink"]], 
                                   with_alleles = FALSE)
  } else {
    r <- ieugwasr::ld_matrix(variant_list, pop = pop, with_alleles = FALSE)
  }
  
  failure <- setdiff(variant_list, colnames(r))
  if (verbose) {
    cat("\nLD information cannot be retrieved for", length(failure), "variants:\n")
    cat(failure, sep="\n")
  }
  
  keep <- intersect(b[["rsid"]], colnames(r))
  ref.keep <- intersect(b[["ref.rsid"]], colnames(r))
  ll <- table(b[["rsid"]], b[["ref.rsid"]])
  ll[,] <- NA
  ll[keep, ref.keep] <- r[keep, ref.keep]
  
  r2 <- sapply(1:nrow(b), function(x) with(b[x, ], ifelse(rsid == ref.rsid, 1, ll[rsid, ref.rsid]^2)))
  
  invisible(cbind(b, r2))
}
