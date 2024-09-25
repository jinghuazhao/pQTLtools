#' Locus novelty check
#'
#' This function checks novelty of a list of loci such as pQTLs against a published list.
#' Both known_loci and query_loci have these variables: chr, pos, uniprot, rsid, prot.
#'
#' @param known_loci A data.frame of published loci.
#' @param query_loci A data.frame of loci whose novelties are unclear.
#' @param ldops arguments for ieugwasr::ld_matrix_local()
#' @param flanking A flanking distance.
#' @param pop The reference population as for ieugwasr::ld_matrix().
#' @param verbose A flag to show nonexistent variants.
#'
#' @return A data.frame containing nonnovel loci.
#' @export
#' @examples
#' \dontrun{
#' options(width=2000)
#' suppressMessages(require(dplyr))
#' suppressMessages(require(openxlsx))
#' # SCALLOP-INF list
#' METAL <- read.delim(file.path(find.package("pQTLtools"),"tests","INF1.METAL")) %>%
#'          dplyr::left_join(gap.datasets::inf1[c("prot","gene")]) %>%
#'          dplyr::mutate(prot=gene,prot_rsid=paste0(uniprot,"-",rsid),chr=Chromosome,pos=Position)
#' # UKB_PPP list
#' require(openxlsx)
#' results <- "/rds/project/jmmh2/rds-jmmh2-results/public/proteomics"
#' url <- file.path(results,"UKB-PPP","doc","sun22.xlsx")
#' ST10 <- read.xlsx(url,"ST10",startRow=3) %>%
#'         dplyr::mutate(uniprot=Target.UniProt,rsid=rsID,prot=Assay.Target) %>%
#'         dplyr::mutate(prot_rsid=paste0(uniprot,"-",rsid))
#' sentinels <- dplyr::left_join(METAL,ST10,by="prot_rsid") %>%
#'              dplyr::select(prot_rsid,cis.trans,rsID) %>%
#'              dplyr::filter(!is.na(rsID))
#' inf1 <- c(with(gap.datasets::inf1,uniprot),with(METAL,uniprot)) %>%
#'         unique()
#' overlap <- dplyr::filter(ST10,uniprot %in% inf1)
#' dim(overlap)
#' UKB_PPP <- dplyr::mutate(overlap,
#'            chrpos=strsplit(overlap[["Variant.ID.(CHROM:GENPOS.(hg37):A0:A1:imp:v1)"]],":"),
#'            chr=as.integer(unlist(lapply(chrpos,"[[",1))),
#'            pos=as.integer(unlist(lapply(chrpos,"[[",2))),
#'            chrpos=paste(chr,pos,sep=":"))
#' # ieugwasr LD reference which requires an up-to-date registration.
#' suppressMessages(require(GenomicRanges))
#' b <- novelty_check(UKB_PPP,METAL)
#' replication <- dplyr::filter(b,r2>=0.8)
#' INF <- "/rds/project/jmmh2/rds-jmmh2-projects/olink_proteomics/scallop/INF/"
#' # write.table(replication,file=file.path(INF,"work","UKB-PPP.txt"),
#' #             row.names=FALSE,quote=FALSE,sep="\t")
#' replication <- read.delim(file.path(find.package("pQTLtools"),"tests","UKB-PPP.txt")) %>%
#'                dplyr::select(known.seqnames,known.rsid,query.rsid,query.prot)
#' variant_list <- unique(c(dplyr::pull(replication,known.rsid),
#'                          dplyr::pull(replication,query.rsid)))
#' load(file.path(find.package("pQTLtools"),"tests","novel_data.rda"))
#' prot_rsid <- with(novel_data,paste0(prot,"-",rsid))
#' prot_rsid_repl <- with(replication,paste0(query.prot,"-",query.rsid))
#' left <- setdiff(prot_rsid,prot_rsid_repl)
#' # local LD reference panel by chromosome
#' # r2 <- LDlinkR::LDmatrix(variant_list,pop="CEU",token=Sys.getenv("LDLINK_TOKEN"))
#' plink <- "/rds/user/jhz22/hpc-work/bin/plink"
#' b <- list()
#' for(i in unique(dplyr::pull(METAL,Chromosome)))
#' {
#'    u <- dplyr::filter(UKB_PPP,chr %in% i) %>%
#'         dplyr::select(chr,pos,uniprot,rsid,prot)
#'    m <- dplyr::filter(METAL,Chromosome %in% i) %>%
#'         dplyr::select(chr,pos,uniprot,rsid,prot)
#'    bfile <- file.path(INF,"INTERVAL","per_chr",paste0("chr",i))
#'    b[[i]] <- novelty_check(u,m,ldops=list(bfile=bfile,plink=plink))
#' }
#' replication2 <- filter(bind_rows(b), r2>=0.8)
#' prot_rsid <- with(novel_data %>%
#'              dplyr::left_join(gap.datasets::inf1[c("prot","gene")]),paste0(gene,"-",rsid))
#' prot_rsid_repl <- with(replication2,paste0(query.prot,"-",query.rsid))
#' novel <- setdiff(prot_rsid,prot_rsid_repl)
#' }

novelty_check <- function(known_loci,query_loci,ldops=NULL,flanking=1e6,pop="EUR",verbose=TRUE)
{
  rsid <- seqnames <- start <- strand <- width <- NA
  bfile <- plink <- NA
  query <- with(known_loci,GenomicRanges::GRanges(seqnames=chr,IRanges::IRanges(start=pos,width=1),
                                                  uniprot=uniprot,rsid=rsid,pos=pos,prot=prot))
  subject <- with(query_loci,GenomicRanges::GRanges(seqnames=chr,IRanges::IRanges(start=pos-flanking,end=pos+flanking),
                                                    uniprot=uniprot,rsid=rsid,pos=pos,prot=prot))
  fo <- GenomicRanges::findOverlaps(query,subject) %>%
        data.frame()
  eq <- subject[fo$subjectHits,]$uniprot==query[fo$queryHits,]$uniprot
  ov <- data.frame(fo[eq,])
  ov1 <- query[ov$queryHits,] %>%
         data.frame() %>%
         select(-strand,-width) %>%
         mutate(rsid=if_else(rsid=="-",paste0("chr",seqnames,":",start),rsid))
  ov2 <- subject[ov$subjectHits,] %>%
         data.frame() %>%
         select(-strand,-width)
  b <- bind_cols(data.frame(ov1) %>% setNames(paste("known",names(ov1),sep=".")),
                 data.frame(ov2) %>% setNames(paste("query",names(ov2),sep=".")))
  variant_list <- unique(c(b[["known.rsid"]],b[["query.rsid"]]))
  if (!is.null(ldops))
  {
    r <- ieugwasr::ld_matrix_local(variants=c(b[["known.rsid"]],b[["query.rsid"]]),
                                   bfile=ldops[["bfile"]],
                                   plink_bin=ldops[["plink"]],
                                   with_alleles=FALSE)
  } else r <- ieugwasr::ld_matrix(variant_list,pop=pop,with_alleles=FALSE)
  failure <- setdiff(variant_list,colnames(r))
  if (verbose) {
     cat("\nLD information cannot be retrieved for", length(failure), "variants:\n")
     cat(failure,sep="\n")
  }
  known.keep <- intersect(b[["known.rsid"]],colnames(r))
  query.keep <- intersect(b[["query.rsid"]],colnames(r))
  ll <- table(b[["known.rsid"]],b[["query.rsid"]])
  ll[,] <- NA
  ll[known.keep,query.keep] <- r[known.keep,query.keep]
  r2 <- sapply(1:nrow(b), function(x) with(b[x, ], ifelse(known.rsid == query.rsid, 1, ll[known.rsid, query.rsid]^2)))
  invisible(mutate(b,r2=r2))
}

# l <- matrix(NA,length(b[["known.rsid"]]),length(b[["query.rsid"]]),dimnames=list(b[["known.rsid"]], b[["query.rsid"]]))
# l[known.keep,query.keep] <- r[known.keep,query.keep]
# r2 <-  sapply(1:nrow(b),function(x) with(b[x,],ifelse(known.rsid==query.rsid,1,l[known.rsid,query.rsid]^2)))

# wget https://www.biorxiv.org/content/biorxiv/early/2022/06/18/2022.06.17.496443/DC2/embed/media-2.xlsx -O sun22.xlsx

