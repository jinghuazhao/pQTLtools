#' Import OpenGWAS
#'
#' A function which imports OpenGWAS.
#'
#' @param opengwas_id An OpenGWAS id.
#' @param region chr:start-end.
#' @param method Method to extract GWAS data.
#' @param verbose Extra information.
#' @param ... Parameters to pass to TwoSampleMR outcome extraction.
#'
#' @details
#' By default, method="TwoSampleMR" should work in all cases with some controls over variant filtering. If a VCF file, \insertCite{lyon21;textual}{pQTLtools}, is known to exist, one can specify method="gwasvcf" to extract a chunk of data.
#'
#' @export
#' @return
#' A summary statistic object. With method="TwoSampleMR" the result is in TwoSampleMR outcome format.
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso \code{\link[TwoSampleMR]{extract_outcome_data}}
#' @examples
#' \dontrun{
#' options(width=200)
#' # method="TwoSampleMR"
#' # GSMR data preparation for Crohn's disease in the LTA region
#' opengwas_id <- "ebi-a-GCST004132"
#' region <- "6:30539831-32542101"
#' n <- 2/(1/12194 + 1/28072)
#' od <- pQTLtools::import_OpenGWAS(opengwas_id,region) %>%
#'       dplyr::distinct() %>%
#'       dplyr::mutate(snpid=gap::chr_pos_a1_a2(chr,pos,effect_allele.outcome,other_allele.outcome),
#'                     effect_allele.outcome=toupper(effect_allele.outcome),
#'                     other_allele.outcome=toupper(other_allele.outcome)) %>%
#'       dplyr::select(snpid,effect_allele.outcome,other_allele.outcome,eaf.outcome,
#'                     beta.outcome,se.outcome,pval.outcome,samplesize.outcome) %>%
#'       setNames(c("SNP","A1","A2","freq","b","se","p","N")) %>%
#'       dplyr::group_by(SNP) %>%
#'       dplyr::slice(which.min(p)) %>%
#'       data.frame()
#' unlink("ebi-a-GCST007432.vcf.gz.tbi")
#' od[is.na(od$N),"N"] <- n
#' write.table(od,quote=FALSE,row.names=FALSE)
#' # method="gwasvcf"
#' gwasvcf::set_bcftools(path=file.path(HPC_WORK,"bin","bcftools"))
#' # MPV ARHGEF3 region
#' opengwas_id <- "ebi-a-GCST004599"
#' region <- "3:56649749-57049749"
#' mpv_ARHGEF3 <- import_OpenGWAS(opengwas_id,region,method="gwasvcf")
#' # all immune-related
#' INF <- Sys.getenv("INF")
#' HPC_WORK <- Sys.getenv("HPC_WORK")
#' opengwas_ids <- scan(file.path(INF,"OpenGWAS","ieu.list"),what="")
#' unavail <- c("ieu-b-18","finn-a-M13_SLE","finn-a-D3_SARCOIDOSIS")
#' opengwas_ids <- subset(opengwas_ids,!opengwas_ids %in% unavail)
#' region <- "1:100-2000000"
#' summary_list = purrr::map(opengwas_ids[1:2], ~import_OpenGWAS(., region, method="gwasvcf"))
#' }
#' @note
#' Adapted function.
#' @keywords utilities

import_OpenGWAS <- function(opengwas_id, region, method="TwoSampleMR", verbose = TRUE, ...)
{
  if (method=="TwoSampleMR") TwoSampleMR::extract_outcome_data(region,opengwas_id,...) else if (method=="gwasvcf")
  {
    opengwas_root <- "https://gwas.mrcieu.ac.uk/files"
    file_path <- paste(opengwas_root,opengwas_id,paste0(opengwas_id,".vcf.gz"),sep="/")
    if(verbose) print(file_path)
  # pending on its exclusive/web use later
  # fetch_table = seqminer::tabix.read.table(tabixFile = file_path, tabixRange = region, stringsAsFactors = FALSE)
  # only possible locally but its GRanges conversion is helpful
  # gwas_stats <- gwasvcf::query_gwas(basename(file_path),chrompos = region)
  # gwas_sumstats <- gwasvcf::vcf_to_granges(gwas_stats)
  # One that does work at the moment
    VariantAnnotation::VcfFile(file_path)
    VariantAnnotation::vcfFields(file_path)
    chr_start_end <- unlist(strsplit(gsub(":", "-", region),"-"))
    seqnames <- chr_start_end[1]
    start <- as.integer(chr_start_end[2])
    end <- as.integer(chr_start_end[3])
    rngs <- GenomicRanges::GRanges(seqnames, IRanges::IRanges(start,end))
    param <- VariantAnnotation::ScanVcfParam(which=rngs)
    vcf <- VariantAnnotation::readVcf(file_path, "hg19", param)
    gwasvcf::vcf_to_granges(vcf)
  }
}
