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
#' INF1_merge <- merge(inf1,
#'                   read.delim(file.path(INF,"work","INF1.merge-rsid"),as.is=TRUE),
#'                   by="prot")
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
#' options(width=250)
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
#' query <- noquote(inf1[["uniprot"]])
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

#' Import eQTL Catalogue
#'
#' An adopted function which imports eQTL Catalogue.
#'
#' @param ftp_path URL.
#' @param region chr:start-end.
#' @param selected_gene_id An Ensembl gene ID.
#' @param column_names Column names of the dataset.
#' @param verbose Extra information.
#'
#' @details
#' This function is based on the eQTL-Catalogue-resources, \insertCite{kerimov21;textual}{pQTLtools}.
#'
#' @export
#' @return A summary statistic object.
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' \dontrun{
#' library(pQTLtools)
#' suppressMessages(invisible(lapply(c("dplyr", "ggplot2", "readr", "coloc",
#'                                     "GenomicRanges","seqminer"),
#'                  require, character.only = TRUE)))
#' # Largely deprecated so b./c. below are local version
#' fp <- file.path(path.package("pQTLtools"),"eQTL-Catalogue","tabix_ftp_paths.tsv")
#' tabix_paths <- read.delim(fp, stringsAsFactors = FALSE) %>% dplyr::as_tibble()
#' # MPV association at the ARHGEF3 locus
#' region <- "3:56615721-57015721"
#' ensGene <- "ENSG00000163947"
#' platelet_df <- dplyr::filter(tabix_paths, study == "CEDAR", tissue_label == "platelet")
#' hdr <- file.path(path.package("pQTLtools"),"eQTL-Catalogue","column_names.CEDAR")
#' column_names <- names(read.delim(hdr))
#' summary_stats <- import_eQTLCatalogue(platelet_df$ftp_path, region,
#'                                       selected_gene_id = ensGene, column_names)
#' summary_stats
#' ggplot(summary_stats, aes(x = position, y = -log(pvalue, 10))) + geom_point()
#' # gwasvcf::set_bcftools(path=file.path(HPC_WORK,"bin","bcftools"))
#' # GWAS sumstat from the same region
#' # manually download and parse with gwasvcf
#' # wget https://gwas.mrcieu.ac.uk/files/ebi-a-GCST004599/ebi-a-GCST004599.vcf.gz
#' # wget https://gwas.mrcieu.ac.uk/files/ebi-a-GCST004599/ebi-a-GCST004599.vcf.gz.tbi
#' # gwas_stats <- gwasvcf::query_gwas("ebi-a-GCST004599.vcf.gz", chrompos = "3:56649749-57049749")
#' # gwas_stats <- gwasvcf::vcf_to_granges(gwas_stats) %>%
#'                 keepSeqlevels("3") %>%
#'                 renameSeqlevels("chr3")
#' # via import_OpenGWAS
#' opengwas_id <- "ebi-a-GCST004599"
#' region <- "3:56649749-57049749"
#' gwas_stats <- import_OpenGWAS(opengwas_id,region) %>%
#'               keepSeqlevels("3") %>%
#'               renameSeqlevels("chr3")
#' f <- file.path(path.package("pQTLtools"),"eQTL-Catalogue","hg19ToHg38.over.chain")
#' chain <- rtracklayer::import.chain(f)
#' gwas_stats_hg38 <- rtracklayer::liftOver(gwas_stats, chain) %>%
#'   unlist() %>%
#'   dplyr::as_tibble() %>%
#'   dplyr::transmute(chromosome = seqnames,
#'                    position = start, REF, ALT, AF, ES, SE, LP, SS) %>%
#'   dplyr::mutate(id = paste(chromosome, position, sep = ":")) %>%
#'   dplyr::mutate(MAF = pmin(AF, 1-AF)) %>%
#'   dplyr::group_by(id) %>%
#'   dplyr::mutate(row_count = n()) %>%
#'   dplyr::ungroup() %>%
#'   dplyr::filter(row_count == 1) %>%
#'   dplyr::mutate(chromosome=gsub("chr","",chromosome))
#' ggplot(gwas_stats_hg38, aes(x = position, y = LP)) + geom_point()
#' # Colocalisation
#' res <- run_coloc(summary_stats, gwas_stats_hg38)
#'
#' # a. all other eQTL datasets
#' microarray_df <- dplyr::filter(tabix_paths, quant_method == "microarray") %>%
#'                  dplyr::mutate(qtl_id = paste(study, qtl_group, sep = "_"))
#' ftp_path_list <- setNames(as.list(microarray_df$ftp_path), microarray_df$qtl_id[1])
#' hdr <- file.path(path.package("pQTLtools"),"eQTL-Catalogue","column_names.CEDAR")
#' column_names <- names(read.delim(hdr))
#' summary_list <- purrr::map(ftp_path_list, ~import_eQTLCatalogue(., region,
#'                            selected_gene_id = ensGene, column_names))
#' coloc_df_microarray <- purrr::map_df(summary_list, ~run_coloc(., gwas_stats_hg38),
#'                                      .id = "qtl_id")
#'
#' # b. Uniformly processed RNA-seq datasets
#' # rnaseq_df <- dplyr::filter(tabix_paths, quant_method == "ge") %>%
#' #              dplyr::mutate(qtl_id = paste(study, qtl_group, sep = "_"))
#' # ftp_path_list <- setNames(as.list(rnaseq_df$ftp_path), rnaseq_df$qtl_id)
#' # hdr <- file.path(path.package("pQTLtools"),"eQTL-Catalogue","column_names.Alasoo")
#' fp <- file.path(find.package("pQTLtools"),"eQTL-Catalogue",
#'                 "tabix_ftp_paths_ge.tsv")
#' # eQTL Catalogue site (deprecated)
#' imported_tabix_paths <- read.delim(fp, stringsAsFactors = FALSE) %>%
#'                         dplyr::as_tibble()
#' # local downloads
#' imported_tabix_paths <- within(imported_tabix_paths,
#'       {
#'         f <- lapply(strsplit(ftp_path,"/csv/|/ge/"),"[",3)
#'         ftp_path <- paste0("~/rds/public_databases/eQTLCatalogue/",f)
#'       })
#' ftp_path_list <- setNames(as.list(imported_tabix_paths$ftp_path),
#'                           imported_tabix_paths$unique_id)
#' hdr <- file.path(path.package("pQTLtools"),"eQTL-Catalogue",
#'                  "column_names.Alasoo")
#' column_names <- names(read.delim(hdr))
#' safe_import <- purrr::safely(import_eQTLCatalogue)
#' summary_list <- purrr::map(ftp_path_list, ~safe_import(., region,
#'                            selected_gene_id = ensGene, column_names))
#' result_list <- purrr::map(summary_list, ~.$result)
#' result_list <- result_list[!unlist(purrr::map(result_list, is.null))]
#' coloc_df_rnaseq <- purrr::map_df(result_list, ~run_coloc(., gwas_stats_hg38),
#'                                  .id = "qtl_id")
#'
#' # c. GTEx_v8 imported eQTL datasets
#' # rnaseq_df <- dplyr::filter(imported_tabix_paths, quant_method == "ge") %>%
#' #              dplyr::mutate(qtl_id = paste(study, qtl_group, sep = "_"))
#' # ftp_path_list <- setNames(as.list(rnaseq_df$ftp_path), rnaseq_df$qtl_id)
#' fp <- file.path(find.package("pQTLtools"),"eQTL-Catalogue",
#'                 "tabix_ftp_paths_gtex.tsv")
#' # eQTL Catalogue site
#' imported_tabix_paths <- read.delim(fp, stringsAsFactors = FALSE) %>%
#'                         dplyr::as_tibble()
#' # local downloads
#' imported_tabix_paths <- within(imported_tabix_paths,
#'        {
#'          f <- lapply(strsplit(ftp_path,"/csv/|/ge/"),"[",3);
#'          ftp_path <- paste0("~/rds/public_databases/GTEx/csv/",f)
#'        })
#' gtex_df <- dplyr::filter(imported_tabix_paths, quant_method == "ge") %>%
#'            dplyr::mutate(qtl_id = paste(study, qtl_group, sep = "_"))
#' ftp_path_list <- setNames(as.list(gtex_df$ftp_path), gtex_df$qtl_id)
#' hdr <- file.path(path.package("pQTLtools"),"eQTL-Catalogue","column_names.GTEx")
#' column_names <- names(read.delim(hdr))
#' safe_import <- purrr::safely(import_eQTLCatalogue)
#' summary_list <- purrr::map(ftp_path_list, ~safe_import(., region,
#'                            selected_gene_id = ensGene, column_names))
#' result_list <- purrr::map(summary_list, ~.$result)
#' result_list <- result_list[!unlist(purrr::map(result_list, is.null))]
#' result_filtered <- purrr::map(result_list, ~dplyr::filter(., !is.na(se)))
#' coloc_df_imported <- purrr::map_df(result_filtered, ~run_coloc(., gwas_stats_hg38),
#'                                    .id = "qtl_id")
#'
#' coloc_df = dplyr::bind_rows(coloc_df_microarray, coloc_df_rnaseq, coloc_df_imported)
#' dplyr::arrange(coloc_df, -PP.H4.abf)
#' ggplot(coloc_df, aes(x = PP.H4.abf)) + geom_histogram()
#' }
#' @note
#' Adapted function.
#' @keywords utilities

import_eQTLCatalogue <- function(ftp_path, region, selected_gene_id, column_names, verbose = TRUE)
{
  if(verbose) print(ftp_path)
  gene_id <- rsid <- chromosome <- position <- id <- n <- row_count <- NULL
  fetch_table <- seqminer::tabix.read.table(tabixFile = ftp_path, tabixRange = region, stringsAsFactors = FALSE) %>% dplyr::as_tibble()
  colnames(fetch_table) <- column_names
  summary_stats <- dplyr::filter(fetch_table, gene_id == selected_gene_id) %>%
    dplyr::select(-rsid) %>%
    dplyr::distinct() %>%
    dplyr::mutate(id = paste(chromosome, position, sep = ":")) %>%
    dplyr::group_by(id) %>%
    dplyr::mutate(row_count = n()) %>% dplyr::ungroup() %>%
    dplyr::filter(row_count == 1)
}

#' Colocalisation analysis
#'
#' The function takes eQTL and GWAS summary statistics for a colocalisation analysis.
#'
#' @param eqtl_sumstats eQTL summary data.
#' @param gwas_sumstats GWAS summary data.
#' @param harmonise a flag to harmonise data.
#' @md
#' @export
#' @return Summary from `coloc.abf`.

run_coloc <- function(eqtl_sumstats, gwas_sumstats, harmonise=TRUE)
{
  if (harmonise)
  {
    chromosome <- position <- ref <- alt <- REF <- ALT <- beta <- se <- ES <- SE <- NA
    eqtl_sumstats <- dplyr::filter(eqtl_sumstats, !any(is.na(c(chromosome,position,ref,alt)))) %>%
                     dplyr::mutate(snpid = gap::chr_pos_a1_a2(chromosome, position, ref, alt))
    gwas_sumstats <- dplyr::filter(gwas_sumstats, !any(is.na(c(chromosome,position,REF,ALT)))) %>%
                     dplyr::mutate(snpid = gap::chr_pos_a1_a2(chromosome, position, REF, ALT)) %>%
                     dplyr::select(-chromosome,-position)
    harmonise_data <- dplyr::left_join(eqtl_sumstats,gwas_sumstats,by='snpid') %>%
                      dplyr::mutate(flag=if_else(ref==REF,1,-1)) %>%
                      dplyr::filter(!is.na(beta)&!is.na(se)&!is.na(ES)&!is.na(SE))
    eqtl_dataset <- with(harmonise_data, list(beta = flag*beta, varbeta = se^2, N = an/2, MAF = maf, snp = snpid, type = "quant"))
    gwas_dataset <- with(harmonise_data, list(beta = ES, varbeta = SE^2, MAF = MAF, N = SS, snp = snpid, type = "quant"))
  } else {
    eqtl_dataset <- list(beta = eqtl_sumstats$beta,
                         varbeta = eqtl_sumstats$se^2,
                         N = (eqtl_sumstats$an)[1]/2,
                         MAF = eqtl_sumstats$maf,
                         type = "quant",
                         snp = eqtl_sumstats$id)
    gwas_dataset <- list(beta = gwas_sumstats$ES,
                         varbeta = gwas_sumstats$SE^2,
                         type = "quant",
                         snp = gwas_sumstats$id,
                         MAF = gwas_sumstats$MAF,
                         N = gwas_sumstats$SS)
  }
  for (p in "coloc") {
    if (length(grep(paste("^package:", p, "$", sep=""), search())) == 0) {
       if (!requireNamespace(p, quietly = TRUE))
           warning(paste("This function needs package `", p, "' to be fully functional; please install", sep=""))
    }
  }
  coloc_res <- coloc::coloc.abf(dataset1 = eqtl_dataset, dataset2 = gwas_dataset, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
  res_formatted <- dplyr::as_tibble(t(as.data.frame(coloc_res$summary)))
}

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
#' library(pQTLtools)
#' # method="TwoSampleMR"
#' # GSMR data preparation for Crohn's disease in the LTA region
#' opengwas_id <- "ebi-a-GCST004132"
#' region <- "6:30539831-32542101"
#' n <- 2/(1/12194 + 1/28072)
#' require(dplyr)
#' od <- import_OpenGWAS(opengwas_id,region) %>%
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

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param dat PARAM_DESCRIPTION, Default: NULL
#' @param type PARAM_DESCRIPTION, Default: 'exposure'
#' @param snps PARAM_DESCRIPTION, Default: NULL
#' @param header PARAM_DESCRIPTION, Default: TRUE
#' @param phenotype_col PARAM_DESCRIPTION, Default: 'Phenotype'
#' @param snp_col PARAM_DESCRIPTION, Default: 'SNP'
#' @param beta_col PARAM_DESCRIPTION, Default: 'beta'
#' @param se_col PARAM_DESCRIPTION, Default: 'se'
#' @param eaf_col PARAM_DESCRIPTION, Default: 'eaf'
#' @param effect_allele_col PARAM_DESCRIPTION, Default: 'effect_allele'
#' @param other_allele_col PARAM_DESCRIPTION, Default: 'other_allele'
#' @param pval_col PARAM_DESCRIPTION, Default: 'pval'
#' @param units_col PARAM_DESCRIPTION, Default: 'units'
#' @param ncase_col PARAM_DESCRIPTION, Default: 'ncase'
#' @param ncontrol_col PARAM_DESCRIPTION, Default: 'ncontrol'
#' @param samplesize_col PARAM_DESCRIPTION, Default: 'samplesize'
#' @param gene_col PARAM_DESCRIPTION, Default: 'gene'
#' @param id_col PARAM_DESCRIPTION, Default: 'id'
#' @param min_pval PARAM_DESCRIPTION, Default: 1e-200
#' @param z_col PARAM_DESCRIPTION, Default: 'z'
#' @param info_col PARAM_DESCRIPTION, Default: 'info'
#' @param chr_col PARAM_DESCRIPTION, Default: 'chr'
#' @param pos_col PARAM_DESCRIPTION, Default: 'pos'
#' @param log_pval PARAM_DESCRIPTION, Default: FALSE
#' @return OUTPUT_DESCRIPTION
#' @rdname format_data.args
#' @keywords internal

format_data.args <- function(dat=NULL, type = "exposure", snps = NULL, header = TRUE,
                             phenotype_col = "Phenotype", snp_col = "SNP", beta_col = "beta",
                             se_col = "se", eaf_col = "eaf", effect_allele_col = "effect_allele",
                             other_allele_col = "other_allele", pval_col = "pval", units_col = "units",
                             ncase_col = "ncase", ncontrol_col = "ncontrol", samplesize_col = "samplesize",
                             gene_col = "gene", id_col = "id", min_pval = 1e-200, z_col = "z",
                             info_col = "info", chr_col = "chr", pos_col = "pos", log_pval = FALSE)
  invisible(list(dat = dat, type = type, snps = snps, header = header,
                 phenotype_col = phenotype_col, snp_col = snp_col, beta_col = beta_col,
                 se_col = se_col, eaf_col = eaf_col, effect_allele_col = effect_allele_col,
                 other_allele_col = other_allele_col, pval_col = pval_col, units_col = units_col,
                 ncase_col = ncase_col, ncontrol_col = ncontrol_col, samplesize_col = samplesize_col,
                 gene_col = gene_col, id_col = id_col, min_pval = min_pval, z_col = z_col,
                 info_col = info_col, chr_col = chr_col, pos_col = pos_col, log_pval = FALSE))

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param snps PARAM_DESCRIPTION, Default: NULL
#' @param outcomes PARAM_DESCRIPTION, Default: NULL
#' @param proxies PARAM_DESCRIPTION, Default: TRUE
#' @param rsq PARAM_DESCRIPTION, Default: 0.8
#' @param align_alleles PARAM_DESCRIPTION, Default: 1
#' @param palindromes PARAM_DESCRIPTION, Default: 1
#' @param maf_threshold PARAM_DESCRIPTION, Default: 0.3
#' @param access_token PARAM_DESCRIPTION, Default: ieugwasr::check_access_token()
#' @param splitsize PARAM_DESCRIPTION, Default: 10000
#' @param proxy_splitsize PARAM_DESCRIPTION, Default: 500
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @seealso
#'  \code{\link[ieugwasr]{check_access_token}}
#' @rdname extract_outcome_data.args
#' @keywords internal

extract_outcome_data.args <- function(snps=NULL, outcomes=NULL, proxies=TRUE, rsq=0.8, align_alleles=1,
                                      palindromes=1, maf_threshold=0.3, access_token=ieugwasr::check_access_token(),
                                      splitsize=10000, proxy_splitsize=500)
  invisible(list(snps=snps,outcomes=outcomes,proxies=proxies,rsq=rsq,align_alleles=align_alleles,
                           palindromes=palindromes,maf_threshold=maf_threshold,access_token=access_token,
                           splitsize=splitsize,proxy_splitsize=proxy_splitsize))

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param dat PARAM_DESCRIPTION, Default: NULL
#' @param clump_kb PARAM_DESCRIPTION, Default: 10000
#' @param clump_r2 PARAM_DESCRIPTION, Default: 0.001
#' @param clump_p1 PARAM_DESCRIPTION, Default: 1
#' @param clump_p2 PARAM_DESCRIPTION, Default: 1
#' @param pop PARAM_DESCRIPTION, Default: 'EUR'
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @rdname clump_data.args
#' @keywords internal

clump_data.args <- function(dat=NULL, clump_kb = 10000, clump_r2 = 0.001, clump_p1 = 1, clump_p2 = 1, pop = "EUR")
  invisible(list(dat=dat,clump_kb=clump_kb,clump_r2=clump_r2,clump_p1=clump_p1,clump_p2=clump_p2,pop=pop))

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param exposure_dat PARAM_DESCRIPTION, Default: NULL
#' @param outcome_dat PARAM_DESCRIPTION, Default: NULL
#' @param action PARAM_DESCRIPTION, Default: 2
#' @return OUTPUT_DESCRIPTION
#' @rdname harmonise_data.args
#' @keywords internal

harmonise_data.args <- function(exposure_dat=NULL, outcome_dat=NULL, action = 2)
  invisible(list(exporesure_dat=exposure_dat,outcome_dat=outcome_dat,action=action))

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param x PARAM_DESCRIPTION
#' @param y PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @keywords internal

swap <- function(x,y)
   eval(parse(text = paste("swap_unique_var_a <-", substitute(x), ";",
   substitute(x), "<-", substitute(y), ";",
   substitute(y), "<-swap_unique_var_a")), envir=parent.frame())

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param message PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @seealso
#'  \code{\link[ggplot2]{ggplot}},\code{\link[ggplot2]{geom_label}},\code{\link[ggplot2]{aes}},\code{\link[ggplot2]{labs}},\code{\link[ggplot2]{theme}},\code{\link[ggplot2]{margin}}
#' @rdname blank_plot
#' @keywords internal

blank_plot <- function(message)
{
   for (p in "ggplot2") {
       if (length(grep(paste("^package:", p, "$", sep=""), search())) == 0) {
          if (!requireNamespace(p, quietly = TRUE))
             warning(paste("This function needs package `", p, "' to be fully functional; please install", sep=""))
       }
   }
   a <- b <- NA
   ggplot2::ggplot(data.frame(a=0,b=0,n=message)) +
   ggplot2::geom_text(ggplot2::aes(x=a,y=b,label=n)) +
   ggplot2::labs(x=NULL,y=NULL) +
   ggplot2::theme(axis.text=ggplot2::element_blank(), axis.ticks=ggplot2::element_blank())
}

#' @title MR scatter plot
#' @description Adaptation of TwoSampleMR::mr_scatter_plot()
#' @param mr_results MR results.
#' @param dat data.
#' @param alpha Confidentce level, Default: 0.05.
#' @return A graphic object.
#' @seealso
#'  \code{\link[plyr]{dlply}},\code{\link[plyr]{mutate}}
#'  \code{\link[TwoSampleMR]{mr_egger_regression}},\code{\link[TwoSampleMR]{mr_egger_regression_bootstrap}}
#'  \code{\link[ggplot2]{ggplot}},\code{\link[ggplot2]{aes}},\code{\link[ggplot2]{geom_crossbar}},\code{\link[ggplot2]{geom_errorbarh}},\code{\link[ggplot2]{geom_point}},\code{\link[ggplot2]{ggtheme}},\code{\link[ggplot2]{geom_abline}},\code{\link[ggplot2]{scale_manual}},\code{\link[ggplot2]{labs}},\code{\link[ggplot2]{theme}},\code{\link[ggplot2]{guide_legend}}
#'  \code{\link[cowplot]{theme_cowplot}}
#' @export
#' @keywords internal

mr_scatter_plot2 <- function (mr_results, dat, alpha=0.05)
{
   for (p in c("ggplot2","plyr")) {
     if (length(grep(paste("^package:", p, "$", sep=""), search())) == 0) {
        if (!requireNamespace(p, quietly = TRUE))
            warning(paste("This function needs package `", p, "' to be fully functional; please install", sep=""))
     }
   }
   c <- qnorm(alpha/2, lower.tail = FALSE)
   a <- b <- id.exposure <- id.outcome <- method <- mr_keep <- se.exposure <- se.outcome <- NA
   mrres <- plyr::dlply(dat, c("id.exposure", "id.outcome"),
      function(d) {
          d <- plyr::mutate(d)
          if (nrow(d) < 2 | with(d,sum(mr_keep)) == 0) return(blank_plot("Insufficient number of SNPs"))
          d <- subset(d, mr_keep)
          index <- with(d, beta.exposure < 0)
          d <- within(d, {beta.exposure[index] <- beta.exposure[index] * -1})
          d <- within(d, {beta.outcome[index] <- beta.outcome[index] * -1})
          mrres <- subset(mr_results, id.exposure == with(d, id.exposure)[1] & id.outcome == with(d, id.outcome)[1])
          mrres <- within(mrres, {a <- 0})
          if ("MR Egger" %in% with(mrres, method)) {
              temp <- with(d, TwoSampleMR::mr_egger_regression(beta.exposure, beta.outcome, se.exposure, se.outcome, default_parameters()))
              mrres <- within(mrres, {a[method == "MR Egger"] <- with(temp, b_i)})
          }
          if ("MR Egger (bootstrap)" %in% with(mrres, method)) {
              temp <- with(d, TwoSampleMR::mr_egger_regression_bootstrap(beta.exposure, beta.outcome, se.exposure, se.outcome, default_parameters()))
              mrres <- within(mrres, {a[method == "MR Egger (bootstrap)"] <- with(temp, b_i)})
          }
          colours <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c",
                       "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928")
          ggplot2::ggplot(data = d, 
          ggplot2::aes(x = beta.exposure, y = beta.outcome)) + 
          ggplot2::geom_errorbar(ggplot2::aes(ymin = beta.outcome - c*se.outcome, ymax = beta.outcome + c*se.outcome), width = 0) +
          ggplot2::geom_errorbarh(ggplot2::aes(xmin = beta.exposure - c*se.exposure, xmax = beta.exposure + c*se.exposure), height = 0) +
          ggplot2::geom_point() +
          ggplot2::theme_bw() +
          cowplot::theme_cowplot(12) + 
          ggplot2::geom_abline(data = mrres, ggplot2::aes(intercept = a, slope = b, colour = method), size = 1, show.legend = TRUE) +
          ggplot2::scale_colour_manual(values = colours) +
          ggplot2::labs(colour = "MR Test", x = paste("SNP effect on", with(d,exposure)[1]), y = paste("SNP effect on", with(d, outcome)[1])) + 
          ggplot2::theme(legend.position = "bottom", legend.direction = "vertical") +
          ggplot2::guides(colour = ggplot2::guide_legend(ncol = 4)) +
          ggplot2::geom_abline(intercept = 0, slope = 0, size = 1) +
          ggplot2::geom_vline(xintercept = 0, size = 1)
      })
   mrres
}

#' @title MR forest plot
#' @description Adaptation from TwoSampleMR::mr_forest_plot()
#' @param singlesnp_results Results based on single variants.
#' @param exponentiate Logic variable to indicate exponentiation, Default: FALSE.
#' @param alpha Confidence level, Default: 0.05.
#' @return A graphic object.
#' @seealso
#'  \code{\link[plyr]{dlply}},\code{\link[plyr]{mutate}}
#'  \code{\link[ggplot2]{ggplot}},\code{\link[ggplot2]{aes}},\code{\link[ggplot2]{ggtheme}},\code{\link[ggplot2]{geom_abline}},\code{\link[ggplot2]{geom_errorbarh}},\code{\link[ggplot2]{geom_point}},\code{\link[ggplot2]{scale_manual}},\code{\link[ggplot2]{theme}},\code{\link[ggplot2]{margin}},\code{\link[ggplot2]{labs}}
#' @export
#' @keywords internal

mr_forest_plot2 <- function (singlesnp_results, exponentiate = FALSE, alpha = 0.05)
{
   for (p in c("ggplot2","plyr")) {
     if (length(grep(paste("^package:", p, "$", sep=""), search())) == 0) {
        if (!requireNamespace(p, quietly = TRUE))
            warning(paste("This function needs package `", p, "' to be fully functional; please install", sep=""))
     }
   }
   c <- qnorm(alpha/2, lower.tail = FALSE)
   b <- id.exposure <- id.outcome <- method <- mr_keep <- se <- se.exposure <- se.outcome <- NA
   res <- plyr::dlply(singlesnp_results, c("id.exposure", "id.outcome"),
          function(d) {
          d <- plyr::mutate(d)
          if (with(d, sum(!grepl("All", SNP))) < 2) return(blank_plot("Insufficient number of SNPs"))
          d <- within(d, {levels(SNP)[levels(SNP) == "All - Inverse variance weighted"] <- "All - IVW"})
          d <- within(d, {levels(SNP)[levels(SNP) == "All - MR Egger"] <- "All - Egger"})
          am <- with(d, grep("All", SNP, value = TRUE))
          d <- within(d, {up <- b + c * se})
          d <- within(d, {lo <- b - c * se})
          d <- within(d, {tot <- 0.01})
          d <- within(d, {tot[SNP %in% am] <- 1})
          d <- within(d, {SNP <- as.character(SNP)})
          nom <- with(d, SNP[!SNP %in% am])
          nom <- with(d, nom[order(b)])
          d <- rbind(d, d[nrow(d), ])
          d <- within(d, {SNP[nrow(d) - 1] <- ""})
          d <- within(d, {b[nrow(d) - 1] <- NA})
          d <- within(d, {up[nrow(d) - 1] <- NA})
          d <- within(d, {lo[nrow(d) - 1] <- NA})
          d <- within(d, {SNP <- ordered(SNP, levels = c(am, "", nom))})
          xint <- 0
          if (exponentiate) {
              d <- within(d, {b <- exp(b)})
              d <- within(d, {up <- exp(up)})
              d <- within(d, {lo <- exp(lo)})
              xint <- 1
          }
          ggplot2::ggplot(d, ggplot2::aes(y = SNP, x = b)) +
          ggplot2::theme_bw() +
#         cowplot::theme_cowplot(12) +
          ggplot2::geom_vline(xintercept = xint, linetype = "dotted") +
          ggplot2::geom_errorbarh(ggplot2::aes(xmin = lo, xmax = up, size = as.factor(tot), colour = as.factor(tot)), height = 0) +
          ggplot2::geom_point(ggplot2::aes(colour = as.factor(tot))) +
          ggplot2::geom_hline(ggplot2::aes(yintercept = which(levels(SNP) %in% ""))) +
          ggplot2::scale_colour_manual(values = c("black", "red")) +
          ggplot2::scale_size_manual(values = c(0.3, 1)) +
          ggplot2::theme(legend.position = "none", axis.text.y = ggplot2::element_text(size = 8),
                         axis.ticks.y = ggplot2::element_line(size = 0), axis.title.x = ggplot2::element_text(size = 8),
                         panel.grid=ggplot2::element_blank()) +
          ggplot2::labs(y = "", x = paste0("MR effect size for\n'", with(d, exposure)[1], "' on '", with(d, outcome)[1], "'"))
          })
   res
}

#' @title MR funnel plot
#' @description Adaptation from TwoSampleMR::mr_funnel_plot()
#' @param singlesnp_results Results based on single variants.
#' @return A graphic object.
#' @seealso
#'  \code{\link[plyr]{dlply}},\code{\link[plyr]{mutate}}
#'  \code{\link[ggplot2]{ggplot}},\code{\link[ggplot2]{aes}},\code{\link[ggplot2]{geom_point}},\code{\link[ggplot2]{ggtheme}},\code{\link[ggplot2]{geom_abline}},\code{\link[ggplot2]{scale_manual}},\code{\link[ggplot2]{labs}},\code{\link[ggplot2]{theme}}
#'  \code{\link[cowplot]{theme_cowplot}}
#' @export
#' @keywords internal

mr_funnel_plot2 <- function (singlesnp_results)
{
  for (p in c("ggplot2","plyr")) {
      if (length(grep(paste("^package:", p, "$", sep=""), search())) == 0) {
         if (!requireNamespace(p, quietly = TRUE))
            warning(paste("This function needs package `", p, "' to be fully functional; please install", sep=""))
      }
  }
  SNP <- b <- se <- NA
  res <- plyr::dlply(singlesnp_results, c("id.exposure", "id.outcome"),
      function(d) {
          d <- plyr::mutate(d)
          if (with(d, sum(!grepl("All", SNP))) < 2) return(blank_plot("Insufficient number of SNPs"))
          am <- with(d, grep("All", SNP, value = TRUE))
          d <- within(d, {SNP <- gsub("All - ", "", SNP)})
          am <- gsub("All - ", "", am)
          ggplot2::ggplot(subset(d, !SNP %in% am), ggplot2::aes(y = 1/se, x = b)) +
          ggplot2::geom_point() +
          ggplot2::theme_bw() +
          cowplot::theme_cowplot(12) +
          ggplot2::geom_vline(data = subset(d, SNP %in% am), ggplot2::aes(xintercept = b, colour = SNP)) +
          ggplot2::scale_colour_manual(values = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c",
                                                  "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928")) +
          ggplot2::labs(y = expression(1/SE[IV]), x = expression(beta[IV]), colour = "MR Method") +
          ggplot2::theme(legend.position = "bottom", legend.direction = "vertical")
      })
  res
}

#' @title MR leave-one-out analysis
#' @description Adapatation from TwoSampleMR::mr_leaveoneout_plot()
#' @param leaveoneout_results Results from leave-one-out analysis.
#' @param alpha PARAM_DESCRIPTION, Default: 0.05.
#' @return A graphic object.
#' @seealso
#'  \code{\link[plyr]{dlply}},\code{\link[plyr]{mutate}}
#'  \code{\link[ggplot2]{ggplot}},\code{\link[ggplot2]{aes}},\code{\link[ggplot2]{ggtheme}},\code{\link[ggplot2]{geom_abline}},\code{\link[ggplot2]{geom_errorbarh}},\code{\link[ggplot2]{geom_point}},\code{\link[ggplot2]{scale_manual}},\code{\link[ggplot2]{theme}},\code{\link[ggplot2]{margin}},\code{\link[ggplot2]{labs}}
#' @export
#' @keywords internal

mr_leaveoneout_plot2 <- function (leaveoneout_results, alpha = 0.05)
{
  for (p in c("ggplot2","plyr")) {
      if (length(grep(paste("^package:", p, "$", sep=""), search())) == 0) {
         if (!requireNamespace(p, quietly = TRUE))
            warning(paste("This function needs package `", p, "' to be fully functional; please install", sep=""))
     }
  }
  c <- qnorm(alpha/2, lower.tail = FALSE)
  a <- b <- id.exposure <- id.outcome <- method <- mr_keep <- se <- se.exposure <- se.outcome <- NA
  res <- plyr::dlply(leaveoneout_results, c("id.exposure", "id.outcome"), function(d) {
      d <- plyr::mutate(d)
      if (with(d, sum(!grepl("All", SNP))) < 3) return(blank_plot("Insufficient number of SNPs"))
      d <- within(d, {up <- b + c * se})
      d <- within(d, {lo <- b - c * se})
      d <- within(d, {tot <- 1})
      d <- within(d, {tot[SNP != "All"] <- 0.01})
      d <- within(d, {SNP <- as.character(SNP)})
      nom <- with(d, SNP[SNP != "All"])
      nom <- nom[order(with(d,b))]
      d <- rbind(d, d[nrow(d), ])
      d <- within(d, {SNP[nrow(d) - 1] <- ""})
      d <- within(d, {b[nrow(d) - 1] <- NA})
      d <- within(d, {up[nrow(d) - 1] <- NA})
      d <- within(d, {lo[nrow(d) - 1] <- NA})
      d <- within(d, {SNP <- ordered(SNP, levels = c("All", "", nom))})
      ggplot2::ggplot(d, ggplot2::aes(y = SNP, x = b)) +
      ggplot2::theme_bw() +
#     cowplot::theme_cowplot(12) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dotted") +
      ggplot2::geom_errorbarh(ggplot2::aes(xmin = lo, xmax = up, size = as.factor(tot), colour = as.factor(tot)), height = 0) +
      ggplot2::geom_point(ggplot2::aes(colour = as.factor(tot))) +
      ggplot2::geom_hline(ggplot2::aes(yintercept = which(levels(SNP) %in% ""))) +
      ggplot2::scale_colour_manual(values = c("black", "red")) +
      ggplot2::scale_size_manual(values = c(0.3, 1)) +
      ggplot2::theme(legend.position = "none", axis.text.y = ggplot2::element_text(size = 8),
                     axis.ticks.y = ggplot2::element_line(size = 0),
                     axis.title.x = ggplot2::element_text(size = 8),
                     panel.grid=ggplot2::element_blank()) +
      ggplot2::labs(y = "", x = paste0("MR leave-one-out sensitivity analysis for\n'", with(d, exposure)[1], "' on '", with(d, outcome)[1], "'"))
  })
  res
}

#' Basic pQTL-MR analysis
#'
#' This function takes exposure and outcome data as produced by `format_data()` and used to perform
#' MR analysis against a list of outcomes; in the latter case it can be data from MR-Base, e.g.,
#' `outcome <- extract_outcome_data(snps=with(exposure,SNP),outcomes=c("ieu-a-7","ebi-a-GCST007432"))`.
#'
#' @param exposure exposure data.
#' @param outcome the counterpart for outcome.
#' @param mr_plot to produce plots.
#' @param prefix a prefix for output files.
#' @param reverse if TRUE, perform reverse MR.
#'
#' @details
#' This function is adapted from \insertCite{zheng20;textual}{pQTLtools} and in
#' spirit similar to `run_TwoSampleMR`.
#' @md
#'
#' @export
#' @return
#' No value is returned but several files.
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' library(TwoSampleMR)
#' library(pQTLtools)
#' fi <- file.path(find.package("pQTLtools",lib.loc=.libPaths()),"tests","Ins.csv")
#' exposure <- format_data(read.csv(fi),type="exposure")
#' fo <- file.path(find.package("pQTLtools",lib.loc=.libPaths()),"tests","Out.csv")
#' outcome <- format_data(read.csv(fo),type="outcome")
#' pqtlMR(exposure, outcome, prefix="IL6R-")
#' pqtlMR(exposure, outcome, prefix="IL6R_rev-",reverse=TRUE)
#' unlink(c("IL6R*","pQTL-combined*"))
#' # Phenotype,SNP,effect_allele,other_allele,eaf,beta,se,pval
#' # ABO,rs505922,C,T,0.313,1.298,0.014,1.2e-1828
#' # LIFR,rs635634,T,C,0.180,-0.300,0.032,6.00E-21
#' # f <- file.path(find.package("pQTLtools",lib.loc=.libPaths()),"tests","ms.ins")
#' # exposure <- format_data(read.table(f, header=TRUE), samplesize_col="N")
#' # SNP Phenotype effect_allele other_allele eaf beta se pval N
#' # rs1800693 TNFB T C 0.6033 0.0282  0.0136 0.0389045   11787
#' # rs2364485 TNFB A C 0.1645 6514963 0.1759 1.62181e-20 11344
#' # https://raw.githubusercontent.com/MRCIEU/epigraphdb-pqtl/master/scripts/MR-pQTL-script.R
#'
#' @note
#' Adapted from script by Jie Zheng.
#' @keywords utilities

pqtlMR <- function(exposure, outcome, mr_plot=FALSE, prefix="pQTL-combined-", reverse=FALSE)
{
   for (p in c("TwoSampleMR")) {
     if (length(grep(paste("^package:", p, "$", sep=""), search())) == 0) {
        if (!requireNamespace(p, quietly = TRUE))
            warning(paste("This function needs package `", p, "' to be fully functional; please install", sep=""))
     }
   }
   id.exposure <- NA
   id.outcome <- NA
   effect_allele.exposure <- NA
   effect_allele.outcome <- NA
   other_allele.exposure <- NA
   other_allele.outcome <- NA
   eaf.exposure <- NA
   eaf.outcome <- NA
   samplesize.exposure <- NA
   beta.exposure <- NA
   beta.outcome <- NA
   se.exposure <- NA
   se.outcome <- NA
   pval.exposure <- NA
   pval.outcome <- NA
   swap_unique_var_a <- NA
   harmonise <- TwoSampleMR::harmonise_data(exposure,outcome)
   if (reverse) harmonise <- subset(within(harmonise,
   {
      swap(id.exposure,id.outcome)
      swap(effect_allele.exposure,effect_allele.outcome)
      swap(other_allele.exposure,other_allele.outcome)
      swap(eaf.exposure,eaf.outcome)
      if(exists("samplesize.exposure") & exists("samplesize.outcome")) swap(samplesize.exposure,samplesize.outcome)
      if(!exists("samplesize.exposure") & exists("samplesize.outcome")) samplesize.outcome <- NA
      swap(beta.exposure,beta.outcome)
      swap(se.exposure,se.outcome)
      swap(pval.exposure,pval.outcome)
    }),select=-swap_unique_var_a)
   result <- single <- NULL
   try(result <- TwoSampleMR::mr(harmonise,method_list=c("mr_wald_ratio", "mr_ivw")))
   try(single <- TwoSampleMR::mr_singlesnp(harmonise))
   invisible(lapply(c("harmonise","result","single"), function(x) if (exists(x))
             write.table(format(get(x),digits=3),file=paste0(prefix,x,".txt"),
                         quote=FALSE,row.names=FALSE,sep="\t")))
   if (mr_plot)
   {
      scatter <- TwoSampleMR::mr_scatter_plot(result, harmonise)
      forest <- TwoSampleMR::mr_forest_plot(single)
      funnel <- TwoSampleMR::mr_funnel_plot(single)
      invisible(sapply(c(scatter,forest,funnel), function(x) print(x)))
   }
}

#' Basic TwoSampleMR analysis
#'
#' Given harmonised data, this function conducts a two-sample MR analysis.
#'
#' @param TwoSampleMRinput Harmonised data.
#' @param mr_plot one of "None", "TwoSampleMR", "pQTLtools" for no, the original and the revised plots, respectively.
#' @param prefix a prefix for output files.
#'
#' @details
#' As TwoSampleMR faces seemingly perplexing options, this function intends to simplify various steps in a two-sample MR
#' as in \insertCite{dt18;textual}{pQTLtools}. It is particularly useful when a large numbher of MRs are necessary,
#' e.g., multiple proteins and their cis/trans regions need to be examined, in which case prefix could direct the output
#' to various directories.
#'
#' @export
#' @return
#' No value is returned but several files.
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' library(TwoSampleMR)
#' library(pQTLtools)
#' outcomes <- "ebi-a-GCST007432"
#' prot <- "MMP.10"
#' type <- "cis"
#' f <- paste0(prot,"-",type,".mrx")
#' d <- read.table(file.path(find.package("pQTLtools",lib.loc=.libPaths()),"tests",f),
#'                 header=TRUE)
#' exposure <- TwoSampleMR::format_data(within(d,{P=10^logP}), phenotype_col="prot", snp_col="rsid",
#'                                      chr_col="Chromosome", pos_col="Posistion",
#'                                      effect_allele_col="Allele1", other_allele_col="Allele2",
#'                                      eaf_col="Freq1", beta_col="Effect", se_col="StdErr",
#'                                      pval_col="P", log_pval=FALSE,
#'                                      samplesize_col="N")
#' clump <- exposure[sample(1:nrow(exposure),nrow(exposure)/80),] # TwoSampleMR::clump_data(exposure)
#' outcome <- TwoSampleMR::extract_outcome_data(snps=exposure$SNP,outcomes=outcomes)
#' harmonise <- TwoSampleMR::harmonise_data(clump,outcome)
#' prefix <- paste(outcomes,prot,type,sep="-")
#' run_TwoSampleMR(harmonise, mr_plot="pQTLtools", prefix=prefix)
#' caption <- "Table. MMP.10 variants and FEV1"
#' knitr::kable(read.delim(paste0(prefix,"-result.txt"),header=TRUE),
#'              caption=paste(caption, "(result)"))
#' knitr::kable(read.delim(paste0(prefix,"-heterogeneity.txt"),header=TRUE),
#'              caption=paste(caption,"(heterogeneity)"))
#' knitr::kable(read.delim(paste0(prefix,"-pleiotropy.txt"),header=TRUE),
#'              caption=paste(caption,"(pleiotropy)"))
#' knitr::kable(read.delim(paste0(prefix,"-single.txt"),header=TRUE),
#'              caption=paste(caption,"(single)"))
#' knitr::kable(read.delim(paste0(prefix,"-loo.txt"),header=TRUE),
#'              caption=paste(caption,"(loo)"))
#' for (x in c("result","heterogeneity","pleiotropy","single","loo"))
#'     unlink(paste0(prefix,"-",x,".txt"))
#'
#' @keywords utilities

run_TwoSampleMR <- function(TwoSampleMRinput, mr_plot="None", prefix="")
{
  for (p in "TwoSampleMR") {
    if (length(grep(paste("^package:", p, "$", sep=""), search())) == 0) {
       if (!requireNamespace(p, quietly = TRUE))
           warning(paste("This function needs package `", p, "' to be fully functional; please install", sep=""))
    }
  }
  harmonise <- TwoSampleMRinput
  result <- heterogeneity <- pleiotropy <- single <- loo <- NULL
  result <- TwoSampleMR::mr(harmonise)
  heterogeneity <- TwoSampleMR::mr_heterogeneity(harmonise)
  pleiotropy <- TwoSampleMR::mr_pleiotropy_test(harmonise)
  single <- TwoSampleMR::mr_singlesnp(harmonise)
  loo <- TwoSampleMR::mr_leaveoneout(harmonise)
  invisible(lapply(c("result","heterogeneity","pleiotropy","single","loo"), function(x)
                   if (exists(x)) write.table(format(get(x),digits=3),
                                              file=paste0(prefix,"-",x,".txt"),
                                              quote=FALSE,row.names=FALSE,sep="\t")))
  type <- match.arg(mr_plot, c("None", "TwoSampleMR", "pQTLtools"))
  if (type == 1) return(0) else if (type == "TwoSampleMR")
  { 
    scatter = TwoSampleMR::mr_scatter_plot(result, harmonise)
    forest <- TwoSampleMR::mr_forest_plot(single)
    funnel <- TwoSampleMR::mr_funnel_plot(single)
    leaveoneout <- TwoSampleMR::mr_leaveoneout_plot(loo)
  } else {
    scatter = mr_scatter_plot2(result, harmonise)
    forest <- mr_forest_plot2(single)
    funnel <- mr_funnel_plot2(single)
    leaveoneout <- mr_leaveoneout_plot2(loo)
  }
  invisible(sapply(c(scatter,forest,funnel,leaveoneout), function(x) print(x)))
}

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

#' Locus novelty check
#'
#' This function checks novelty of a list of loci such as pQTLs against a published list.
#' Both known_loci and query_loci have these variables: chr, pos, uniprot, rsid, prot.
#'
#' @param known_loci A data.frame of published loci.
#' @param query_loci A data.frame of loci whose novelties are unclear.
#' @param flanking A flanking distance.
#' @param pop The reference population as for ieugwasr::ld_matrix().
#' @param verbose A flag to show nonexistent variants.
#'
#' @return A data.frame containing nonnovel loci.
#' @export
#' @examples
#' \dontrun{ 
#'  options(width=2000)
#'  suppressMessages(require(dplyr))
#'  # SCALLOP-INF list
#'  METAL <- read.delim("~/pQTLtools/tests/INF1.METAL") %>%
#'           mutate(prot_rsid=paste0(uniprot,"-",rsid),pos=Position)
#'  # UKB_PPP list
#'  require(openxlsx)
#'  url <- "~/rds/results/public/proteomics/UKB-PPP/sun22.xlsx"
#'  ST10 <- read.xlsx(url,"ST10",startRow=3) %>%
#'          mutate(uniprot=Target.UniProt,rsid=rsID,prot=Assay.Target) %>%
#'          mutate(prot_rsid=paste0(uniprot,"-",rsid))
#'  sentinels <- left_join(METAL,ST10,by="prot_rsid") %>%
#'               select(prot_rsid,cis.trans,rsID) %>%
#'               filter(!is.na(rsID))
#'  inf1 <- c(with(pQTLtools::inf1,uniprot),with(METAL,uniprot)) %>%
#'          unique()
#'  overlap <- filter(ST10,uniprot %in% inf1)
#'  dim(overlap)
#'  UKB_PPP <- mutate(overlap,
#'             chrpos=strsplit(overlap[["Variant.ID.(CHROM:GENPOS.(hg37):A0:A1:imp:v1)"]],":"),
#'             chr=lapply(chrpos,"[[",1),
#'             pos=lapply(chrpos,"[[",2),
#'             chrpos=paste(chr,pos,sep=":"))
#'  suppressMessages(require(GenomicRanges))
#'  b <- novelty_check(UKB_PPP,METAL)
#'  replication <- filter(b,r2>=0.8)
#'  # Mutual uses of pQTLtools/SCALLOP-INF
#'  write.table(replication,file="~/INF/work/UKB-PPP.txt",row.names=FALSE,quote=FALSE,sep="\t")
#'  load("~/pQTLtools/tests/novel_data.rda")
#'  prot_rsid <- with(novel_data,paste0(prot,"-",rsid))
#'  prot_rsid_repl <- with(replication,paste0(query.prot,"-",query.rsid))
#'  left <- setdiff(prot_rsid,prot_rsid_repl)
#' }

novelty_check <- function(known_loci,query_loci,flanking=1e6,pop="EUR",verbose=TRUE)
{
  rsid <- seqnames <- start <- strand <- width <- NA
  query <- with(known_loci,GenomicRanges::GRanges(seqnames=as.integer(chr),IRanges::IRanges(start=as.integer(pos),width=1),
                                                  uniprot=uniprot,rsid=rsid,prot=prot))
  subject <- with(query_loci,GenomicRanges::GRanges(seqnames=Chromosome,IRanges::IRanges(start=pos-flanking,end=pos+flanking),
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
  b <- bind_cols(data.frame(ov1),data.frame(ov2))
  names(b) <- c(paste("known",names(ov1),sep="."),paste("query",names(ov2),sep="."))
  variant_list <- unique(c(b[["known.rsid"]],b[["query.rsid"]]))
  r <- ieugwasr::ld_matrix(variant_list,pop=pop,with_alleles=FALSE)
  failure <- setdiff(variant_list,colnames(r))
  if (verbose) {
     cat("\nLD information cannot be retrieved for", length(failure), "variants:\n")
     cat(failure,sep="\n")
  }
  known.keep <- intersect(b[["known.rsid"]],colnames(r))
  query.keep <- intersect(b[["query.rsid"]],colnames(r))
  ll <- table(b$known.rsid,b$query.rsid)
  ll[,] <- NA
  ll[known.keep,query.keep] <- r[known.keep,query.keep]
  r2 <- sapply(1:nrow(b), function(x) with(b[x, ], ifelse(known.rsid == query.rsid, 1, ll[known.rsid, query.rsid]^2)))
  invisible(mutate(b,r2=r2))
}

# l <- matrix(NA,length(b[["known.rsid"]]),length(b[["query.rsid"]]),dimnames=list(b[["known.rsid"]], b[["query.rsid"]]))
# l[known.keep,query.keep] <- r[known.keep,query.keep]
# r2 <-  sapply(1:nrow(b),function(x) with(b[x,],ifelse(known.rsid==query.rsid,1,l[known.rsid,query.rsid]^2)))

# wget https://www.biorxiv.org/content/biorxiv/early/2022/06/18/2022.06.17.496443/DC2/embed/media-2.xlsx -O sun22.xlsx

#' QTL lookup
#'
#' This function takes MR results (involving pQTL/QTL) to look up QTLs in trait GWASs given a P value and
#' linkage disequilibrium (LD) cutoffs.
#' The rationale is that a pQTL may not necessarily be in strong LD with QTL but some other independent signal
#' in the same region. The interest then lands on association signals below a given P value (p_threshold) and above
#' a given LD (r2_threshold) thresholds.
#'
#' @param d directory where `dat` (below) is held.
#' @param dat MR results (`protein`, `id`, `pqtl`,`p`, `qtl`, `p_qtl`) whose `proxy`, `p_proxy` and `rsq` variables will be updated.
#' @param panel reference panel.
#' @param p_threshold cutoff of QTL association from trait GWASs.
#' @param r2_threshold cutoff of r^2.
#' @param pop reference population from 1000Genomes for LD calculation.
#' @param plink_bin PLINK executable file whose binary files are indicated in `bfile` variable of `dat`.
#' @param r when specified, the LD(r) is output.
#' @param r2 when specified, the LD(r^2) is output.
#' @param xlsx a non-null specification indicates name of an output Excel workbook.
#' 
#' @return A data.frame containing the looked up loci.
#' @export
#' @examples
#' \dontrun{
#' options(width=200)
#' INF <- Sys.getenv("INF")
#' suppressMessages(library(dplyr))
#' d <- file.path(INF,"mr","gsmr","trait")
#' inf1 <- select(gap.datasets::inf1,prot,target.short)
#' gsmr_efo <- read.delim(file.path(INF,"mr","gsmr","gsmr-efo.txt")) %>%
#'             left_join(inf1,by=c('protein'='target.short')) %>%
#'             mutate(file_gwas=paste(prot,id,"rsid.txt",sep="-"),
#'                    bfile=file.path(INF,"INTERVAL","per_chr",
#'                                    paste0("interval.imputed.olink.chr_",chr)),
#'                    proxy=NA,p_proxy=NA,rsq=NA)
#' proxies <- qtl_lookup(gsmr_efo,plink_bin="/rds/user/jhz22/hpc-work/bin/plink",
#'                       xlsx=file.path(INF,"mr","gsmr","r2_INTERVAL.xlsx")) %>%
#'            select(protein,id,Disease,fdr,pqtl,p,qtl,p_qtl,proxy,p_proxy,rsq)
#' write.table(proxies,file=file.path(INF,"mr","gsmr","r2_INTERVAL.tsv"),
#'             row.names=FALSE,quote=FALSE,sep="\t")
#' }

qtl_lookup <- function(d,dat,panel="1000Genomes",p_threshold=1e-3,r2_threshold=0.8,pop="EUR",
                       plink_bin=NULL,r=NULL,r2=NULL,xlsx=NULL)
{
  protein <- p <- SNP <- prot <- Disease <- protein <- fdr <- qtl <- p_qtl <- rsq <- NULL
  for(i in 1:nrow(dat))
  {
     z <- dplyr::slice(dat,i)
     pqtl <- z[["pqtl"]]
     cat(z[["file_gwas"]],"\n")
     gwas <- read.table(file.path(d,basename(z[["file_gwas"]])),header=TRUE) %>%
             dplyr::arrange(p)
     h <-  dplyr::filter(gwas,p<=p_threshold)
     panel_snps <- c(pqtl,dplyr::pull(h,SNP))
     if (panel=="1000Genomes")
     {
        if (length(panel_snps)>500) stop("too many SNPs -- put a more stringent ptreshold")
        xx <- ieugwasr::ld_matrix(panel_snps,pop="EUR")
     } else
     xx <- ieugwasr::ld_matrix(panel_snps,with_alleles=TRUE,pop=pop,bfile=z[["bfile"]],plink_bin=plink_bin)
     cn <- colnames(xx)
     inside <- pqtl==gsub("_[A-Z]*","",cn)
     nn <- c(cn[inside],cn[!inside])
     if (length(r)==1) r_mat <- xx else r_mat <- xx[nn,nn]
     if(!is.null(r)) write.table(r_mat,file=paste(prot,id,Disease,pqtl,"r.txt",sep="-"))
     r2_mat <- r_mat^2
     if(!is.null(r2)) write.table(r2_mat,file=paste(prot,id,Disease,pqtl,"r2.txt",sep="-"))
     colnames(r2_mat) <- gsub("_[A-Z]*","",colnames(r2_mat))
     rownames(r2_mat) <- gsub("_[A-Z]*","",rownames(r2_mat))
     snps <- intersect(pull(h,SNP),colnames(r2_mat))
     success <- 1
     if (length(snps)<1) next else if (length(snps)==1) dat[i,c("proxy","p_proxy","rsq")] <- c(z[c("qtl","p_qtl")],1)
     else while(length(snps)>1)
     {
       success <- -1
       proxy <- snps[1]
       snps <- setdiff(snps,proxy)
       r2_i <- r2_mat[z[["pqtl"]],proxy]
       p_proxy <- filter(gwas,SNP==proxy) %>%
                  slice(which.min(p)) %>%
                  pull(p)
       cat("Same locus",i,z[["protein"]],z[["id"]],z[["Disease"]],z[["pqtl"]],z[["qtl"]],proxy,r2_i,z[["p_qtl"]],"\n",sep="\t")
       if(!is.null(r2_i)&!is.na(r2_i)) if(r2_i>r2_threshold) {success <-1; break}
     }
     if (success==-1)  while(length(snps)>1)
     {
       success <- -1
       proxy <- snps[1]
       snps <- setdiff(snps,proxy)
       r2_i <- r2_mat[z[["pqtl"]],proxy]
       p_proxy <- filter(gwas,SNP==proxy) %>%
                  slice(which.min(p)) %>%
                  pull(p)
       cat("Independent locus",i,z[["protein"]],z[["id"]],z[["Disease"]],z[["pqtl"]],z[["qtl"]],proxy,r2_i,z[["p_qtl"]],"\n",sep="\t")
       if(!is.null(r2_i)&!is.na(r2_i)) if(r2_i<r2_threshold) {success <-1; break}
     }
     if (success==1)
     {
       dat[i,"proxy"] <- proxy
       dat[i,"p_proxy"] <- p_proxy
       dat[i,"rsq"] <- r2_i
     }
  }
  if (!is.null(xlsx))
  {
    wb <- openxlsx::createWorkbook(xlsx)
    hs <- openxlsx::createStyle(textDecoration="BOLD", fontColour="#FFFFFF", fontSize=12, fontName="Arial Narrow", fgFill="#4F80BD")
    proxies <- dplyr::select(dat,protein,id,Disease,fdr,pqtl,p,qtl,p_qtl,proxy,p_proxy,rsq)
    for (sheet in "proxies")
    {
      openxlsx::addWorksheet(wb,sheet,zoom=150)
      openxlsx::writeData(wb,sheet,sheet,xy=c(1,1),headerStyle=openxlsx::createStyle(textDecoration="BOLD",
                          fontColour="#FFFFFF", fontSize=14, fontName="Arial Narrow", fgFill="#4F80BD"))
      body <- get(sheet)
      openxlsx::writeDataTable(wb, sheet, body, xy=c(1,2), headerStyle=hs, firstColumn=TRUE, tableStyle="TableStyleMedium2")
      openxlsx::freezePane(wb, sheet, firstActiveCol=3, firstActiveRow=3)
      width_vec <- apply(body, 2, function(x) max(nchar(as.character(x))+2, na.rm=TRUE))
    # width_vec_header <- nchar(colnames(body))+2
      openxlsx::setColWidths(wb, sheet, cols = 1:ncol(body), widths = width_vec)
      openxlsx::writeData(wb, sheet, tail(body,1), xy=c(1, nrow(body)+2), colNames=FALSE, borders="rows", borderStyle="thick")
    }
    openxlsx::saveWorkbook(wb, file=xlsx, overwrite=TRUE)
  }
  dat
}
