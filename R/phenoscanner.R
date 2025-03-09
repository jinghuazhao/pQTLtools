#' phenoscanner genequeries in batches
#'
#' R/phenoscanner only allows for certain number of items supplied. This simple function return
#' a large number of calls in batches as well as generating SNPIDs.
#'
#' @param genelist a list of SNPs.
#' @param catalogue "None","eQTL","mQTL","methQTL","pQTL","GWAS".
#' @param proxies "None", "AFR","AMR","EAS","EUR","SAS".
#' @param p p value threshold.
#' @param r2 r2 for LD.
#' @param build 37, 38.
#' @param wait a flag to wait for 1hr for every 50 genes.
#'
#' @details
#' Batches are generated and queries are combined into one.
#'
#' @export
#' @return The returned value is a list containing genes and results.
#' @seealso \code{\link[phenoscanner]{phenoscanner}}
#' @examples
#' \dontrun{
#' # single gene
#'   genequeries("TNFRSF11B")
#' }
#' @keywords utilities

genequeries <- function(genelist,catalogue="pQTL",proxies="EUR",p=5e-8,r2=0.8,build=37,wait=TRUE)
{
  a1 <- a2 <- hg19_coordinates <- hg38_coordinates <- NULL
  batches <- split(genelist,ceiling(seq_along(genelist)/10))
  g <- r <- vector('list',length(batches))
  for(i in 1:length(batches))
  {
    cat("Block",i,"\n")
    if (wait) if (i%%6==0) Sys.sleep(60*60)
    q <- phenoscanner::phenoscanner(genequery=batches[[i]], catalogue=catalogue,
                                    proxies=proxies, pvalue=p, r2=r2, build=build)
    g[[i]] <- with(q,genes)
    r[[i]] <- with(q,results)
  }
  genes <- do.call("rbind",g)
  results <- within(do.call("rbind",r),
  {
    swap_a1 <- as.character(a1)
    swap_a2 <- as.character(a2)
    swap <- swap_a1 > swap_a2
    a1[swap] <- swap_a2[swap]
    a2[swap] <- swap_a1[swap]
    if (build==37) snpid <- paste0(hg19_coordinates,"_",a1,"_",a2)
    else if (build==38) snpid <- paste0(hg38_coordinates,"_",a1,"_",a2)
    rm(swap_a1,swap_a2,swap)
  })
  list(genes=genes,results=results)
}

#' phenoscanner regionqueries in batches
#'
#' R/phenoscanner only allows for certain number of items supplied. This simple function return
#' a large number of calls in batches as well as generating SNPIDs.
#'
#' @param regionlist a list of SNPs
#' @param catalogue "None","eQTL","mQTL","methQTL","pQTL","GWAS".
#' @param proxies "None", "AFR","AMR","EAS","EUR","SAS".
#' @param p p value threshold.
#' @param r2 r2 for LD.
#' @param build 37, 38.
#' @param wait a flag to wait for 1hr for every 50 regions.
#'
#' @details
#' Batches are generated and queries are combined into one.
#'
#' @export
#' @return The returned value is a list containing tiles, regions and results.
#'
#' @seealso \code{\link[phenoscanner]{phenoscanner}}
#'
#' @examples
#' \dontrun{
#' # single region
#' regionqueries("chr17:26691290-26700110")
#'
#' # SCALLOP -- SomaLogic lookup from PhenoScanner
#' INF <- Sys.getenv("INF")
#' INF1_merge <- merge(gap.datasets::inf1,
#'                     read.delim(file.path(INF,"work","INF1.merge-rsid"),as.is=TRUE),
#'                     by="prot")
#' INF1_merge_uniprot <- with(INF1_merge,unique(uniprot))
#' SomaLogic_INF1_merge <- subset(SomaLogic160410,UniProt %in% INF1_merge_uniprot)
#' regions <- subset(INF1_merge,uniprot %in% with(SomaLogic_INF1_merge,UniProt))
#' singletons <- with(regions, Start-End<=2)
#' flank <- 5e+2
#' regions[singletons,"Start"] <- regions[singletons,"Start"] - flank
#' regions[singletons,"End"] <- regions[singletons,"End"] + flank
#' reset <- with(regions,Start < 0)
#' regions[reset,"Start"] <- 0
#' r <- regionqueries(with(regions,paste0(Chrom,":",Start,"-",End)))
#' save(r,file="INF1_merge.rda",compress='xz')
#' r2 <- with(r,
#' {
#'  region_ext <- cbind(tiles,regions)
#'  results_ext <- merge(region_ext,results,by="region")
#'  ord <- with(results_ext,order(group))
#'  results_ext[ord,]
#' })
#' results <- subset(r2,pmid=="29875488")
#' grp <- names(table(with(results,group)))
#' sink("INF1_merge.txt")
#' for(g in as.numeric(grp))
#' {
#'   uniprot <- regions[g,"uniprot"]
#'   SNP <- regions[g,"SNP"]
#'   print(regions[g,])
#'   s <- subset(results,group==g&rsid==SNP)
#'   vars <- c("region","group","rsid","hg19_coordinates","hgnc","beta","se","p","snpid") 
#'   if(nrow(s)>1) print(s[vars])
#' }
#' sink()
#' }
#' @note
#' adapted from custom codings
#' @keywords utilities

regionqueries <- function(regionlist,catalogue="pQTL",proxies="EUR",p=5e-8,r2=0.8,build=37,wait=TRUE)
{
  a1 <- a2 <- hg19_coordinates <- hg38_coordinates <- NULL
  lrl <- strsplit(regionlist,":|-")
  chr <- as.character(lapply(lrl,"[[",1))
  start <- as.integer(lapply(lrl,"[[",2))
  end <- as.integer(lapply(lrl,"[[",3))
  gr <- GenomicRanges::GRanges(seqnames=chr,IRanges::IRanges(start,end))
  print(GenomicRanges::width(gr))
  tiles <- GenomicRanges::tile(gr, width=1e+6)
  print(GenomicRanges::width(tiles))
  regionlist_ext <- with(as.data.frame(tiles),paste0(seqnames,":",start,"-",end))
  cat("Conducting queries for",length(regionlist_ext),"regions.\n")
  batches <- split(regionlist_ext,ceiling(seq_along(regionlist_ext)/10))
  s <- r <- vector('list',length(batches))
  for(i in 1:length(batches))
  {
    cat("Block",i,"\n")
    if (wait) if (i%%6==0) Sys.sleep(60*60)
    q <- phenoscanner::phenoscanner(regionquery=batches[[i]], catalogue=catalogue,
                                    proxies=proxies, pvalue=p, r2=r2, build=build)
    s[[i]] <- with(q,regions)
    r[[i]] <- with(q,results)
  }
  regions <- do.call("rbind",s)
  results <- within(do.call("rbind",r),
  {
    swap_a1 <- as.character(a1)
    swap_a2 <- as.character(a2)
    swap <- swap_a1 > swap_a2
    a1[swap] <- swap_a2[swap]
    a2[swap] <- swap_a1[swap]
    if (build==37) snpid <- paste0(hg19_coordinates,"_",a1,"_",a2)
    else if (build==38) snpid <- paste0(hg38_coordinates,"_",a1,"_",a2)
    rm(swap_a1,swap_a2,swap)
  })
  list(tiles=tiles,regions=regions,results=results)
}

#' phenoscanner snpqueries in batches
#'
#' R/phenoscanner only allows for certain number of items supplied. This simple function return
#' a large number of calls in batches as well as generating SNPIDs.
#'
#' @param snplist a list of SNPs.
#' @param block_size size of each query block.
#' @param waiting_time time (in seconds) to wait between query blocks.
#' @param catalogue "None","eQTL","mQTL","methQTL","pQTL","GWAS".
#' @param proxies "None", "AFR","AMR","EAS","EUR","SAS".
#' @param p p value threshold.
#' @param r2 r2 for LD.
#' @param build 37, 38.
#' @param wait a flag to wait for 1hr for every 500 SNPs.
#'
#' @details
#' Batches are generated and queries are combined into one.
#'
#' @export
#' @return
#' The returned value is a list containing snps and results:
#' @seealso \code{\link[phenoscanner]{phenoscanner}}
#'
#' @examples
#' \dontrun{
#'  # single SNP
#'  snpqueries("rs704")
#'  # SCALLOP/INF
#' INF <- Sys.getenv("INF")
#' rsid <- scan(paste(INF,'work','INF1.merge.snp',sep='/'),"")
#' r <- snpqueries(rsid,catalogue='pQTL',p=1e-11)
#' INTERVAL_Olink <- subset(with(r,results),efo=='EFO_0004747' & pmid=='29875488')
#' save(INTERVAL_Olink,file='INTERVAL_Olink.rda',compress='xz')
#' # --- query intersect proteins ---
#' # SomaLogic intersect
#' SomaLogic_overlap_list <- subset(st4,UniProt %in% intersect_list)
#' r <- snpqueries(SomaLogic_overlap_list[,6],catalogue='pQTL',p=1e-11)
#' SomaLogic_overlap <- subset(with(r,results),efo=='EFO_0004747' & pmid=='29875488')
#' save(SomaLogic_overlap_list,SomaLogic_overlap,file='SomaLogic_overlap.rda',compress='xz')
#' SomaLogic_result <- merge(SomaLogic_overlap_list,SomaLogic_overlap,
#'                           by.x="Sentinel.variant*",by.y="snp")
#' # Olink intersect
#' INF1_merge_rsid <- read.delim(paste(INF,"work","INF1.merge-rsid",sep="/"))
#' INF1_merge_rsid_uniprot <- merge(INF1_merge_rsid,inf1,by="prot")
#' Olink_overlap_list <- subset(INF1_merge_rsid_uniprot,uniprot %in% intersect_list)
#' r <- snpqueries(with(Olink_overlap_list,MarkerName),catalogue='pQTL',p=1e-11)
#' Olink_overlap <- subset(with(r,results),efo=='EFO_0004747' & pmid=='29875488')
#' save(Olink_overlap_list,Olink_overlap,file='Olink_overlap.rda',compress='xz')
#' Olink_result <- merge(Olink_overlap_list,Olink_overlap,by.x="MarkerName",by.y="snp")
#' }
#' @note
#' adapted from custom codings
#' @keywords utilities

snpqueries <- function(snplist,block_size=100,waiting_time=60*60,
                       catalogue="pQTL",proxies="EUR",p=5e-8,r2=0.8,build=37,wait=TRUE)
{
  a1 <- a2 <- hg19_coordinates <- hg38_coordinates <- NULL
  batches <- split(snplist,ceiling(seq_along(snplist)/block_size))
  s <- r <- vector('list',length(batches))
  for(i in 1:length(batches))
  {
    cat("Block",i,"\n")
    if (wait) if (i%%6==0) Sys.sleep(waiting_time)
    q <- phenoscanner::phenoscanner(snpquery=batches[[i]], catalogue=catalogue,
                                    proxies=proxies, pvalue=p, r2=r2, build=build)
    s[[i]] <- with(q,snps)
    r[[i]] <- with(q,results)
  }
  snps <- within(do.call("rbind",s),
  {
    swap_a1 <- as.character(a1)
    swap_a2 <- as.character(a2)
    swap <- swap_a1 > swap_a2
    a1[swap] <- swap_a2[swap]
    a2[swap] <- swap_a1[swap]
    if (build==37) snpid <- paste0(hg19_coordinates,"_",a1,"_",a2)
    else if (build==38) snpid <- paste0(hg38_coordinates,"_",a1,"_",a2)
    rm(swap_a1,swap_a2,swap)
  })
  results <- within(do.call("rbind",r),
  {
    swap_a1 <- as.character(a1)
    swap_a2 <- as.character(a2)
    swap <- swap_a1 > swap_a2
    a1[swap] <- swap_a2[swap]
    a2[swap] <- swap_a1[swap]
    if (build==37) snpid <- paste0(hg19_coordinates,"_",a1,"_",a2)
    else if (build==38) snpid <- paste0(hg38_coordinates,"_",a1,"_",a2)
    rm(swap_a1,swap_a2,swap)
  })
  list(snps=snps,results=results)
}
