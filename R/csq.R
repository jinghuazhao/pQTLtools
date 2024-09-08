#' Variant consequence
#'
#' This function maps consequences (CSQ) of a given list of variants to those in an annotated
#' set based on facilities such as variant effect predictor (VEP). In the case of
#' protein-altering variants (PAVs), many consequences can be involved. The procedure also
#' takes into account linkage disequilibrium (LD) in flanking regions.
#'
#' @param query_loci A data.frame of loci whose consequences are to be obtained.
#' @param annotated_loci A data.frame of annotated loci.
#' @param consequences A set of consequences to check against.
#' @param ldops arguments for ieugwasr::ld_matrix_local()
#' @param flanking A flanking distance.
#' @param pop The reference population as for ieugwasr::ld_matrix().
#' @param verbose A flag to show nonexistent variants.
#'
#' @return A data.frame with an indicator with/without those consequences.
#' @export
#' @examples
#' \dontrun{
#' options(width=2000)
#' suppressMessages(require(dplyr))
#' suppressMessages(require(stringr))
#' # SCALLOP-INF list
#' METAL <- read.delim(file.path(find.package("pQTLtools"),"tests","INF1.METAL")) %>%
#'          dplyr::left_join(gap.datasets::inf1[c("prot","gene")]) %>%
#'          dplyr::mutate(prot=gene,prot_rsid=paste0(uniprot,"-",rsid),
#'                        chr=Chromosome,pos=Position)
#' # VEP output
#' vep <- "/rds/project/rds-zuZwCZMsS0w/Caprion_proteomics/analysis/bgen/vep"
#' consequences <- c("frameshift_variant","missense_variant","splice_acceptor_variant",
#'                   "splice_acceptor_variant","start_lost","stop_gained")
#' pattern <- paste(consequences, collapse = "|")
#' suppressMessages(require(GenomicRanges))
#' INF <- "/rds/project/rds-zuZwCZMsS0w/olink_proteomics/scallop/INF"
#' plink <- "/rds/user/jhz22/hpc-work/bin/plink"
#' # CSQ
#' b <- list()
#' for(i in unique(pull(METAL,Chromosome)))
#' {
#'    m <- dplyr::filter(METAL,Chromosome %in% i) %>%
#'         dplyr::select(chr,pos,MarkerName,prot) %>%
#'         dplyr::mutate(rsid=gsub("chr","",MarkerName)) %>%
#'         dplyr::select(-MarkerName)
#'    u <- read.delim(file.path(vep,paste0("chr",i,".tab.gz"))) %>%
#'         dplyr::select(Chrom,Pos,X.Uploaded_variation,Consequence) %>%
#'         setNames(c("chr","pos","rsid","csq"))
#'    bfile <- file.path(INF,"INTERVAL","per_chr",paste0("snpid",i))
#'    b[[i]] <- csq(m,u,pattern,ldops=list(bfile=bfile,plink=plink))
#'    names(b[[i]]) <- i
#' }
#' b[["23"]] <- mutate(b[["X"]],seqnames="23",ref.seqnames="23")
#' replication <- dplyr::filter(bind_rows(b[-which(names(b)=="X")]),r2>=0.8) %>%
#'                dplyr::rename(gene=prot,ref.gene=ref.prot) %>%
#'                dplyr::mutate(seqnames=as.integer(seqnames),pos=as.integer(pos)) %>%
#'                dplyr::arrange(seqnames,pos) %>%
#'                dplyr::select(-ref.seqnames,-ref.start,-ref.end,-seqnames,-pos)
#' }

csq <- function(query_loci,annotated_loci,pattern,ldops=NULL,flanking=1e6,pop="EUR",verbose=TRUE)
{
  rsid <- seqnames <- start <- strand <- width <- NA
  bfile <- plink <- NA
  subject <- with(query_loci,
             GenomicRanges::GRanges(seqnames=chr,IRanges::IRanges(start=pos-flanking,end=pos+flanking),
                                    rsid=rsid,pos=pos,prot=prot))
  query <- with(filter(annotated_loci,str_detect(csq, pattern)),
           GenomicRanges::GRanges(seqnames=chr,IRanges::IRanges(start=pos-flanking,end=pos+flanking),
                                  rsid=rsid,pos=pos,csq=csq))
  ov <- GenomicRanges::findOverlaps(query,subject) %>%
        data.frame()
  ov1 <- subject[ov$subjectHits,] %>%
         data.frame() %>%
         select(-strand,-width)
  ov2 <- query[ov$queryHits,] %>%
         data.frame() %>%
         select(-strand,-width) %>%
         mutate(rsid=if_else(rsid=="-",paste0("chr",seqnames,":",start),rsid))
  b <- bind_cols(data.frame(ov1),
                 data.frame(ov2) %>% setNames(paste("ref",names(ov2),sep=".")))
  b_granges <- makeGRangesFromDataFrame(b, keep.extra.columns = TRUE)
  variant_list <- unique(c(b[["rsid"]],b[["ref.rsid"]]))
  if (!is.null(ldops))
  {
    r <- ieugwasr::ld_matrix_local(variants = c(b[["rsid"]], b[["ref.rsid"]]),
                                   bfile = ldops[["bfile"]],
                                   plink_bin = ldops[["plink"]], 
                                   with_alleles = FALSE)
  } else r <- ieugwasr::ld_matrix(variant_list, pop = pop, with_alleles = FALSE)
  failure <- setdiff(variant_list,colnames(r))
  if (verbose) {
     cat("\nLD information cannot be retrieved for", length(failure), "variants:\n")
     cat(failure,sep="\n")
  }
  keep <- intersect(b[["rsid"]],colnames(r))
  ref.keep <- intersect(b[["ref.rsid"]],colnames(r))
  ll <- table(b[["rsid"]],b[["ref.rsid"]])
  ll[,] <- NA
  ll[keep,ref.keep] <- r[keep,ref.keep]
  r2 <- sapply(1:nrow(b), function(x) with(b[x, ], ifelse(rsid == ref.rsid, 1, ll[rsid, ref.rsid]^2)))
  invisible(mutate(b,r2=r2))
}

# gunzip -c chr1.tab.gz | sed "/#'/d" | cut -f9 | sort | uniq | awk '{print "#",$0}'
# 3_prime_UTR_variant
# 3_prime_UTR_variant,NMD_transcript_variant
# 5_prime_UTR_variant
# downstream_gene_variant
# frameshift_variant
# frameshift_variant,splice_region_variant
# frameshift_variant,start_lost,start_retained_variant
# frameshift_variant,stop_retained_variant
# inframe_insertion
# inframe_insertion,splice_region_variant
# intergenic_variant
# intron_variant
# intron_variant,NMD_transcript_variant
# intron_variant,non_coding_transcript_variant
# mature_miRNA_variant
# missense_variant
# missense_variant,NMD_transcript_variant
# missense_variant,splice_region_variant
# missense_variant,splice_region_variant,NMD_transcript_variant
# missense_variant,stop_retained_variant
# non_coding_transcript_exon_variant
# splice_acceptor_variant
# splice_acceptor_variant,non_coding_transcript_variant
# splice_donor_5th_base_variant,intron_variant
# splice_donor_5th_base_variant,intron_variant,non_coding_transcript_variant
# splice_donor_region_variant,intron_variant
# splice_donor_region_variant,intron_variant,non_coding_transcript_variant
# splice_donor_region_variant,non_coding_transcript_exon_variant
# splice_donor_variant
# splice_donor_variant,non_coding_transcript_variant
# splice_polypyrimidine_tract_variant,intron_variant
# splice_polypyrimidine_tract_variant,intron_variant,NMD_transcript_variant
# splice_polypyrimidine_tract_variant,intron_variant,non_coding_transcript_variant
# splice_region_variant,3_prime_UTR_variant
# splice_region_variant,5_prime_UTR_variant
# splice_region_variant,intron_variant
# splice_region_variant,intron_variant,NMD_transcript_variant
# splice_region_variant,intron_variant,non_coding_transcript_variant
# splice_region_variant,non_coding_transcript_exon_variant
# splice_region_variant,non_coding_transcript_variant
# splice_region_variant,splice_polypyrimidine_tract_variant,intron_variant
# splice_region_variant,splice_polypyrimidine_tract_variant,intron_variant,NMD_transcript_variant
# splice_region_variant,splice_polypyrimidine_tract_variant,intron_variant,non_coding_transcript_variant
# splice_region_variant,synonymous_variant
# splice_region_variant,synonymous_variant,NMD_transcript_variant
# start_lost
# stop_gained
# stop_gained,frameshift_variant
# stop_gained,splice_region_variant
# stop_lost
# stop_lost,splice_region_variant
# stop_retained_variant
# synonymous_variant
# synonymous_variant,NMD_transcript_variant
# upstream_gene_variant
