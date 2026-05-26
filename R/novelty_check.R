#' Locus novelty check using LD overlap
#'
#' This function assesses whether query loci are novel or match known loci
#' by combining genomic proximity (± flanking window) and linkage disequilibrium (LD).
#'
#' It supports:
#' - 1000 Genomes LD reference via ieugwasr
#' - Local PLINK reference panels via ld_matrix_local()
#'
#' The function:
#' 1. Finds overlapping loci within a flanking distance
#' 2. Matches loci by gene/protein (uniprot)
#' 3. Computes LD (r and r^2) between known and query variants
#' 4. Returns per-pair LD values for downstream novelty/replication assessment
#'
#' @param known_loci Data.frame of known/published loci.
#' Must contain columns: chr, pos, uniprot, rsid, prot.
#'
#' @param query_loci Data.frame of query loci to evaluate.
#' Must contain columns: chr, pos, uniprot, rsid, prot.
#'
#' @param ldops Optional list specifying local LD computation:
#'   \describe{
#'     \item{bfile}{PLINK binary prefix (bed/bim/fam)}
#'     \item{plink}{Path to PLINK executable}
#'   }
#'
#' @param flanking Genomic window (bp) around query loci used for overlap.
#' Default is 1e6 (±1 Mb).
#'
#' @param pop 1000 Genomes population code (e.g., "EUR") used when ldops = NULL.
#'
#' @param verbose Logical; if TRUE prints missing LD variants.
#'
#' @return A data.frame with paired known/query loci and LD r^2 values.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # 1000G mode
#' novelty_check(known_loci, query_loci)
#'
#' # Local PLINK mode
#' novelty_check(
#'   known_loci,
#'   query_loci,
#'   ldops = list(
#'     bfile = "/path/interval.imputed.olink.chr_3",
#'     plink = "/path/plink"
#'   )
#' )
#' }
#'
novelty_check <- function(
    known_loci,
    query_loci,
    ldops = NULL,
    flanking = 1e6,
    pop = "EUR",
    verbose = TRUE
) {
    required_cols <- c("chr","pos","uniprot","rsid","prot")
    if (!all(required_cols %in% names(known_loci))) {
        stop("known_loci missing required columns")
    }
    if (!all(required_cols %in% names(query_loci))) {
        stop("query_loci missing required columns")
    }
    make_gr <- function(df, flank = 0) {
        GenomicRanges::GRanges(
            seqnames = df$chr,
            ranges = IRanges::IRanges(
                start = pmax(1, df$pos - flank),
                end = df$pos + flank
            ),
            uniprot = df$uniprot,
            rsid = df$rsid,
            prot = df$prot,
            pos = df$pos
        )
    }
    known_gr <- make_gr(known_loci, flank = 0)
    query_gr <- make_gr(query_loci, flank = flanking)
    hits <- GenomicRanges::findOverlaps(known_gr, query_gr)
    if (length(hits) == 0) return(data.frame())
    hits_df <- data.frame(
        known_idx = queryHits(hits),
        query_idx = subjectHits(hits)
    )
    hits_df <- hits_df[
        known_gr[hits_df$known_idx]$uniprot ==
        query_gr[hits_df$query_idx]$uniprot,
    ]
    if (nrow(hits_df) == 0) return(data.frame())
    b <- data.frame(
        known.rsid = known_gr[hits_df$known_idx]$rsid,
        query.rsid = query_gr[hits_df$query_idx]$rsid,
        known.prot = known_gr[hits_df$known_idx]$prot,
        query.prot = query_gr[hits_df$query_idx]$prot,
        known.uniprot = known_gr[hits_df$known_idx]$uniprot,
        query.uniprot = query_gr[hits_df$query_idx]$uniprot,
        stringsAsFactors = FALSE
    )
    variant_list <- unique(c(b$known.rsid, b$query.rsid))
    if (is.null(ldops)) {
        r <- ieugwasr::ld_matrix(
            variants = variant_list,
            pop = pop,
            with_alleles = FALSE
        )
    } else {
        r <- ieugwasr::ld_matrix_local(
            variants = variant_list,
            bfile = ldops$bfile,
            plink_bin = ldops$plink,
            with_alleles = FALSE
        )
    }
    if (!is.matrix(r)) {
        stop("LD matrix computation failed")
    }
    colnames(r) <- gsub("_[A-Z]+$", "", colnames(r))
    rownames(r) <- gsub("_[A-Z]+$", "", rownames(r))
    missing <- setdiff(variant_list, colnames(r))
    if (verbose && length(missing) > 0) {
        message("Missing LD variants: ", length(missing))
    }
    r2_vals <- mapply(
        function(kr, qr) {
            if (kr == qr) return(1)
            if (!(kr %in% colnames(r)) || !(qr %in% rownames(r))) return(NA)
            r[kr, qr]^2
        },
        b$known.rsid,
        b$query.rsid
    )
    b$r2 <- r2_vals
    b
}
