#' phenoscanner genequeries in batches
#'
#' R/phenoscanner only allows for certain number of items supplied. This simple function return
#' a large number of calls in batches as well as generating SNPIDs.
#'
#' @md
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
#' @references
#' Sun BB, et al. (2018). Genomic atlas of the human plasma proteome. *Nature* 558: 73-79.
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
#' @md
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
#' @return The returned value is a list containing tiles, regions and results.
#'
#' @references
#' Sun BB, et al. (2018). Genomic atlas of the human plasma proteome. *Nature* 558: 73-79.
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
#' @md
#' @param snplist a list of SNPs.
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
#' @references
#' Sun BB, et al. (2018). Genomic atlas of the human plasma proteome. *Nature* 558: 73-79.
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

snpqueries <- function(snplist,catalogue="pQTL",proxies="EUR",p=5e-8,r2=0.8,build=37,wait=TRUE)
{
  a1 <- a2 <- hg19_coordinates <- hg38_coordinates <- NULL
  batches <- split(snplist,ceiling(seq_along(snplist)/100))
  s <- r <- vector('list',length(batches))
  for(i in 1:length(batches))
  {
    cat("Block",i,"\n")
    if (wait) if (i%%6==0) Sys.sleep(60*60)
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
#'
#' @return A UniProt-ID mapping
#'
#' @references
#' See https://www.uniprot.org/help/api_idmapping
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
#' @md
#' @param ftp_path URL.
#' @param region chr:start-end.
#' @param selected_gene_id An Ensembl gene ID.
#' @param column_names Column names of the dataset.
#' @param verbose Extra information.
#'
#' @details
#' This function is based on the eQTL-Catalogue-resources.
#'
#' @export
#' @return A summary statistic object.
#'
#' @references
#' Kerimov N., et al. (2020). "eQTL Catalogue: a compendium of uniformly processed human gene expression and splicing QTLs",
#'    bioRxiv: 2020.2001.2029.924266, https://www.ebi.ac.uk/eqtl/.
#'
#' @examples
#' \dontrun{
#' library(pQTLtools)
#' invisible(lapply(c("dplyr", "ggplot2", "readr", "coloc", "GenomicRanges","seqminer"),
#'                  require, character.only = TRUE))
#' ftp <- file.path(path.package("pQTLtools"),"eQTL-Catalogue","tabix_ftp_paths.tsv")
#' tabix_paths <- read.delim(ftp, sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>%
#'                dplyr::as_tibble()
#' tfpi <- file.path(find.package("pQTLtools", lib.loc=.libPaths()),"eQTL-Catalogue",
#'                  "tabix_ftp_paths_gtex.tsv")
#' imported_tabix_paths <- read.delim(tfpi, sep = "\t", stringsAsFactors = FALSE) %>%
#'                         dplyr::as_tibble()
#'
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
#'   renameSeqlevels("3") %>%
#'   dplyr::as_tibble() %>%
#'   dplyr::transmute(chromosome = seqnames,
#'                    position = start, REF, ALT, AF, ES, SE, LP, SS) %>%
#'   dplyr::mutate(id = paste(chromosome, position, sep = ":")) %>%
#'   dplyr::mutate(MAF = pmin(AF, 1-AF)) %>%
#'   dplyr::group_by(id) %>%
#'   dplyr::mutate(row_count = n()) %>%
#'   dplyr::ungroup() %>%
#'   dplyr::filter(row_count == 1)
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
#' rnaseq_df <- dplyr::filter(tabix_paths, quant_method == "ge") %>%
#'              dplyr::mutate(qtl_id = paste(study, qtl_group, sep = "_"))
#' ftp_path_list <- setNames(as.list(rnaseq_df$ftp_path), rnaseq_df$qtl_id)
#' hdr <- file.path(path.package("pQTLtools"),"eQTL-Catalogue","column_names.Alasoo")
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
#' rnaseq_df <- dplyr::filter(imported_tabix_paths, quant_method == "ge") %>%
#'              dplyr::mutate(qtl_id = paste(study, qtl_group, sep = "_"))
#' ftp_path_list <- setNames(as.list(rnaseq_df$ftp_path), rnaseq_df$qtl_id)
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
#' The function takes eQTL and GWAS summary statistics for a colocalisation analysis..
#'
#' @md
#' @param eqtl_sumstats eQTL summary data.
#' @param gwas_sumstats GWAS summary data.
#' @param harmonise a flag to harmonise data.
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
#' @md
#' @param opengwas_id An OpenGWAS id.
#' @param region chr:start-end.
#' @param verbose Extra information
#'
#' @details
#' This function is derived from SCALLOP/INF work.
#'
#' @return
#' A summary statistic object
#'
#' @references
#' Lyon M, Andrews SJ, Elsworth B, Gaunt TR, Hemani G, Marcora E. The variant call format provides efficient and robust storage of GWAS summary statistics. bioRxiv 2020.05.29.115824; doi: https://doi.org/10.1101/2020.05.29.115824
#'
#' @examples
#' \dontrun{
#' options(width=200)
#' gwasvcf::set_bcftools(path=file.path(HPC_WORK,"bin","bcftools"))
#' # MPV ARHGEF3 region
#' opengwas_id <- "ebi-a-GCST004599"
#' region <- "3:56649749-57049749"
#' mpv_ARHGEF3 <- import_OpenGWAS(opengwas_id,region)
#' # all immune-related
#' INF <- Sys.getenv("INF")
#' HPC_WORK <- Sys.getenv("HPC_WORK")
#' opengwas_ids <- scan(file.path(INF,"OpenGWAS","ieu.list"),what="")
#' unavail <-c("ieu-b-18","finn-a-M13_SLE","finn-a-D3_SARCOIDOSIS")
#' opengwas_ids <- subset(opengwas_ids,!opengwas_ids %in% unavail)
#' region <- "1:100-2000000"
#' library(pQTLtools)
#' summary_list = purrr::map(opengwas_ids[1:2], ~import_OpenGWAS(., region))
#' }
#' @note
#' Adapted function.
#' @keywords utilities

import_OpenGWAS <- function(opengwas_id, region, verbose = TRUE)
{
  opengwas_root <- "https://gwas.mrcieu.ac.uk/files"
  file_path <- paste(opengwas_root,opengwas_id,paste0(opengwas_id,".vcf.gz"),sep="/")
  if(verbose) print(file_path)
# pending on its exclusive/web use later
# fetch_table = seqminer::tabix.read.table(tabixFile = file_path, tabixRange = region, stringsAsFactors = FALSE)
# only possible locally but its GRanges conversion is helpful
# gwas_stats <- query_gwas(basename(file_path),chrompos = region)
# gwas_sumstats <- vcf_to_granges(gwas_stats)
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
                     axis.ticks.y = ggplot2::element_line(size = 0), axis.title.x = ggplot2::element_text(size = 8),panel.grid=ggplot2::element_blank()) +
      ggplot2::labs(y = "", x = paste0("MR leave-one-out sensitivity analysis for\n'", with(d, exposure)[1], "' on '", with(d, outcome)[1], "'"))
  })
  res
}

#' Basic pQTL-MR analysis
#'
#' This function takes data intrumental variables as produced by format_data() and
#' used to perform MR analysis against a list of outcomes from MR-Base.
#'
#' @md
#' @param ivs Instrumental variables from format_data().
#' @param ids A list of MR-Base IDs.
#' @param mr_plot to produce plots.
#' @param prefix a prefix for output files.
#' @param reverse if TRUE, perform reverse MR.
#'
#' @details
#' This function is based on TwoSampleMR.
#'
#' @export
#' @return
#' No value is returned but several files.
#'
#' @references
#' Zheng J, et al. (2020). Phenome-wide Mendelian randomization mapping the influence of the plasma proteome on complex diseases. *Nature Genetics* 52(10): 1122-1131.
#'
#' @examples
#' library(TwoSampleMR)
#' library(pQTLtools)
#' # Original examples
#' f <- file.path(find.package("pQTLtools",lib.loc=.libPaths()),"tests","Ins.csv")
#' ivs <- format_data(read.csv(f))
#' ids <- c("ieu-a-7","ebi-a-GCST007432")
#' pqtlMR(ivs, ids, mr_plot=FALSE)
#' # A bidirectional analysis
#' f <- file.path(find.package("pQTLtools",lib.loc=.libPaths()),"tests","ms.ins")
#' ivs <- format_data(read.table(f, header=TRUE), samplesize_col="N")
#' ids <- "ieu-b-18"
#' # MR
#' pqtlMR(ivs, ids, prefix="MS-")
#' # reverse MR
#' pqtlMR(ivs, ids, ,prefix="MS_rev-",reverse=TRUE)
#'
#' @note
#' Adapted from script by Jie Zheng.
#' @keywords utilities

pqtlMR <- function(ivs, ids, mr_plot=FALSE, prefix="pQTL-combined-", reverse=FALSE)
{
   for (p in c("TwoSampleMR")) {
     if (length(grep(paste("^package:", p, "$", sep=""), search())) == 0) {
        if (!requireNamespace(p, quietly = TRUE))
            warning(paste("This function needs package `", p, "' to be fully functional; please install", sep=""))
     }
   }
   outcome_data <- TwoSampleMR::extract_outcome_data(snps=with(ivs,SNP),outcomes=ids)
   harmonise <- TwoSampleMR::harmonise_data(ivs,outcome_data)
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
#  main MR analysis
   try(result <- TwoSampleMR::mr(harmonise,method_list=c("mr_wald_ratio", "mr_ivw")))
#  single SNP MR using Wald ratio
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
#' @md
#' @param TwoSampleMRinput Harmonised data.
#' @param mr_plot one of "None", "TwoSampleMR", "pQTLtools" for no, the original and the revised plots, respectively.
#' @param prefix a prefix for output files.
#'
#' @details
#' As TwoSampleMR faces seemingly perplexing options, this function intends to simplify various steps in a two-sample MR. It is
#' particularly useful when a large numbher of MRs are necessary, e.g., multiple proteins and their cis/trans regions need to be examined,
#' in which case prefix could direct the output to various directories.
#'
#' @export
#' @return
#' No value is returned but several files.
#'
#' @references
#' Dimou NL, Tsilidis KK. A Primer in Mendelian Randomization Methodology with a Focus on Utilizing Published Summary Association Data. In
#' Evangelos Evangelou (ed.), Genetic Epidemiology: Methods and Protocols, Methods in Molecular Biology, vol. 1793,
#' https://doi.org/10.1007/978-1-4939-7868-7_13, Springer Science+Business Media, LLC, part of Springer Nature 2018
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
#' clump <- TwoSampleMR::clump_data(exposure)
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
#' @note
#' Adapted from script by Dimou NL, Tsilidis KK.
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
#' @md
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
#' @md
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

#' @title biomaRt
#' @description Curation of biomaRt
#' @docType data
#' @keywords datasets internal
#' @format A data frame with 42198 rows and 11 variables:
#' \describe{
#'   \item{\code{ensembl_gene_id}}{ENSEMBL gene id}
#'   \item{\code{chromosome_name}}{Chromosome name [1-22,X,Y]}
#'   \item{\code{start_position}}{start}
#'   \item{\code{end_position}}{end}
#'   \item{\code{description}}{Description}
#'   \item{\code{hgnc_symbol}}{HGNC symbol}
#'   \item{\code{ensembl_transcript_id}}{ENSEMBL transcript id}
#'   \item{\code{transcription_start_site}}{TSS}
#'   \item{\code{transcript_start}}{Transcript start}
#'   \item{\code{transcript_end}}{Transcript end}
#'   \item{\code{uniprotswissprot}}{UnitProt id}
#' }
#' @details extraction using R.

"biomaRt"

#' @title Caprion panel
#' @description Information based on Caprion pilot studies
#' @docType data
#' @keywords datasets internal
#' @format A data frame with 987 rows and 7 variables:
#' \describe{
#'   \item{\code{Protein}}{Protein name as in UniProt}
#'   \item{\code{Accession}}{UniProt id}
#'   \item{\code{Gene}}{HGNC symbol}
#'   \item{\code{Protein.Description}}{Detailed information on protein}
#'   \item{\code{GO.Cellular.Component}}{GO Ceullular component}
#'   \item{\code{GO.Function}}{GO function}
#'   \item{\code{GO.Process}}{GO process}
#' }
#' @details See the Caprion repository involving its use.

"caprion"

#' @title hg19 information
#' @description protein information
#' @docType data
#' @keywords datasets internal
#' @format A data frame with 62559 rows and 8 variables:
#' \describe{
#'   \item{\code{chr}}{Chromosome [chr1-22,X,Y,...]}
#'   \item{\code{start}}{start}
#'   \item{\code{end}}{end}
#'   \item{\code{width}}{width}
#'   \item{\code{strand}}{strand}
#'   \item{\code{ENSEMBL}}{ENSMEBL gene}
#'   \item{\code{SYMBOL}}{HGNC symbol}
#'   \item{\code{UNIPROT}}{UniProt id}
#' }
#' @details Curation from R

"hg19"

#' @title hg19 Table
#' @description Gene information
#' @docType data
#' @keywords datasets internal
#' @format A data frame with 19872 rows and 12 variables:
#' \describe{
#'   \item{\code{X.chrom}}{Chromosome}
#'   \item{\code{chromStart}}{chromStart}
#'   \item{\code{chromEnd}}{chromEnd}
#'   \item{\code{strand}}{Strand}
#'   \item{\code{acc}}{UniProt id}
#'   \item{\code{uniprotName}}{Protein}
#'   \item{\code{accList}}{List of UniProt ids}
#'   \item{\code{isoIds}}{isoIds}
#'   \item{\code{geneName}}{geneName}
#'   \item{\code{geneSynonyms}}{geneSynonyms}
#'   \item{\code{hgncSym}}{HGNC symbol}
#'   \item{\code{ensGene}}{ENSEMBL gene}
#' }
#' @details Curation from UCSC.

"hg19Tables"

#' @title Olink/INF1 panel
#' @description The panel is based on SCALLOP-INF
#' @docType data
#' @keywords datasets internal
#' @format A data frame with 92 rows and 9 variables:
#' \describe{
#'   \item{\code{uniprot}}{UniProt id}
#'   \item{\code{prot}}{Protein}
#'   \item{\code{target}}{Protein target name}
#'   \item{\code{target.short}}{Protein target short name}
#'   \item{\code{gene}}{HGNC symbol}
#'   \item{\code{chr}}{chromosome[1-13,16-17,19-22]}
#'   \item{\code{start}}{start}
#'   \item{\code{end}}{end}
#'   \item{\code{ensembl_gene_id}}{ENSEMBL gene}
#' }
#' @details Assembled for SCALLOP-INF

"inf1"

#' @title Olink/NGS panel
#' @description Information based on pilot studies
#' @docType data
#' @keywords datasets internal
#' @format A data frame with 1472 rows and 3 variables:
#' \describe{
#'   \item{\code{UniProt}}{UniProt id}
#'   \item{\code{Assay}}{Experimental assay}
#'   \item{\code{Panel}}{Olink panel}
#' }
#' @details Curated from R.

"Olink_NGS"

#' @title Olink/qPCR panels
#' @description Information on all qPCR panels
#' @docType data
#' @keywords datasets internal
#' @format A data frame with 1112 rows and 7 variables:
#' \describe{
#'   \item{\code{UniProt}}{UniProt id}
#'   \item{\code{Panel}}{Panels}
#'   \item{\code{Target}}{Protein}
#'   \item{\code{gene}}{HGNC symbol}
#'   \item{\code{chr}}{Chromosome}
#'   \item{\code{start}}{start}
#'   \item{\code{end}}{end}
#' }
#' @details Curated from Excel.

"Olink_qPCR"

#' @title Somascan panel
#' @description This is based on panel used in Sun et al. (2018) Nature.
#' @docType data
#' @keywords datasets internal
#' @format A data frame with 5178 rows and 10 variables:
#' \describe{
#'   \item{\code{SOMAMER_ID}}{Somamer id}
#'   \item{\code{UniProt}}{UniProt id}
#'   \item{\code{Target}}{Protein target}
#'   \item{\code{TargetFullName}}{Protein target full name}
#'   \item{\code{chr}}{chromosome[1-22,X,Y]}
#'   \item{\code{start}}{start}
#'   \item{\code{end}}{end}
#'   \item{\code{entGene}}{entrez gene}
#'   \item{\code{ensGene}}{ENSEMBL gene}
#'   \item{\code{extGene}}{external gene}
#' }
#' @details from the INTERVAL study.

"SomaLogic160410"

#' @title SomaScan v4.1
#' @description This is also the latest panel
#' @docType data
#' @keywords datasets internal
#' @format A data frame with 7288 rows and 6 variables:
#' \describe{
#'   \item{\code{#}}{A serial number}
#'   \item{\code{SeqID}}{SeqID}
#'   \item{\code{Human.Target.or.Analyte}}{Human target/analyte}
#'   \item{\code{UniProt.ID}}{UniProt id}
#'   \item{\code{GeneID}}{HGNC symbol}
#'   \item{\code{Type}}{"Protein"}
#' }
#' @details obtained directly from SomaLogic.

"SomaScanV4.1"

#' @title SWATH-MS panel
#' @description Curated during INTERVAL pilot study.
#' @docType data
#' @keywords datasets internal
#' @format A data frame with 684 rows and 5 variables:
#' \describe{
#'   \item{\code{Accession}}{UniProt id}
#'   \item{\code{accList}}{List of UniProt ids}
#'   \item{\code{uniprotName}}{Protein}
#'   \item{\code{ensGene}}{ENSEMBL gene}
#'   \item{\code{geneName}}{HGNC symbol}
#' }
#' @details As above.

"swath_ms"

#' @title Supplementary table 4
#' @description Supplementary information for Sun et al. (2018) Nature
#' @docType data
#' @keywords datasets internal
#' @format A data frame with 1980 rows and 31 variables:
#' \describe{
#'   \item{\code{Locus.ID}}{Locus id}
#'   \item{\code{SOMAmer.ID}}{SOMAmer id}
#'   \item{\code{Target}}{Protein}
#'   \item{\code{Target.fullname}}{Protein full name}
#'   \item{\code{UniProt}}{UniProt id}
#'   \item{\code{Sentinel.variant*}}{Sentinel variant}
#'   \item{\code{Chr}}{Chromosome}
#'   \item{\code{Pos}}{Position}
#'   \item{\code{Region.start}}{Region start}
#'   \item{\code{Region.end}}{Region end}
#'   \item{\code{Effect.Allele.(EA)}}{Effect allele}
#'   \item{\code{Other.Allele.(OA)}}{Other allele}
#'   \item{\code{EAF}}{Effect allele frequency}
#'   \item{\code{INFO}}{Information score}
#'   \item{\code{cis/.trans}}{"cis"/"trans"}
#'   \item{\code{Mapped.gene}}{Mapped gene}
#'   \item{\code{No..conditionally.significant.variants}}{Number of conditionally significant variants}
#'   \item{\code{Previously.reported}}{Previously reported}
#'   \item{\code{Replicates?}}{Yes/No}
#'   \item{\code{beta}}{b}
#'   \item{\code{SE}}{SE}
#'   \item{\code{p}}{p value}
#'   \item{\code{beta}}{b}
#'   \item{\code{SE}}{SE}
#'   \item{\code{p}}{p value}
#'   \item{\code{beta}}{b}
#'   \item{\code{SE}}{SE}
#'   \item{\code{p}}{p value}
#'   \item{\code{Uncorrelated.with.PAV.(r2≥0.1)}}{Uncorrelated with PAV}
#'   \item{\code{Significant.after.adjusting.for.PAVs}}{Significant after adjustment for PAVs}
#'   \item{\code{Is.a.cis.eQTL.for.same.gene?}}{Is a cis eQTL for the same gene?}
#' }
#' @details As above.

"st4"

#' @title Supplementary table 6
#' @description Supplementary information for Sun et al. (2018) Nature
#' @docType data
#' @keywords datasets internal
#' @format A data frame with 163 rows and 20 variables:
#' \describe{
#'   \item{\code{Locus.ID}}{Locus id}
#'   \item{\code{Sentinel.variant*}}{Sentinel variant}
#'   \item{\code{Chr}}{Chromosome [1-12,14-20,22]}
#'   \item{\code{Pos}}{Position}
#'   \item{\code{SOMAmer.ID}}{SOMAmer id}
#'   \item{\code{Target}}{Protein}
#'   \item{\code{Target.fullname}}{Protein full name}
#'   \item{\code{UniProt}}{UniProt id}
#'   \item{\code{cis/.trans}}{"cis"/"trans"}
#'   \item{\code{Mapped.gene}}{Mapped gene}
#'   \item{\code{Effect.Allele.(EA)}}{Effect allele}
#'   \item{\code{Other.Allele.(OA)}}{Other allele}
#'   \item{\code{Previously.reported}}{Previously reported (0/1)}
#'   \item{\code{beta}}{b}
#'   \item{\code{SE}}{SE}
#'   \item{\code{p}}{p value}
#'   \item{\code{beta}}{b}
#'   \item{\code{SE}}{SE}
#'   \item{\code{p}}{p value}
#'   \item{\code{Replicates?}}{Yes/No}
#' }
#' @details As above.

"st6"

#' @title Supplementary table 18
#' @description Supplementary information for Sun et al. (2018) Nature
#' @docType data
#' @keywords datasets internal
#' @format A data frame with 3622 rows and 3 variables:
#' \describe{
#'   \item{\code{Number}}{A serial number}
#'   \item{\code{Analyte}}{Name}
#'   \item{\code{UniProt.ID(s)}}{UniProt id}
#' }
#' @details As above.

"st18"
