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

swap <- function(x,y)
   eval(parse(text = paste("swap_unique_var_a <-", substitute(x), ";",
   substitute(x), "<-", substitute(y), ";",
   substitute(y), "<-swap_unique_var_a")), envir=parent.frame())

pqtlMR <- function(Ins,Ids,prefix="INF1",reverse=FALSE)
{
  exposure <- outcome <- id.exposure <- id.outcome <- effect_allele.exposure <- effect_allele.outcome <- NULL
  other_allele.exposure <- other_allele.outcome <- eaf.exposure <- eaf.outcome <- samplesize.exposure <- NULL
  beta.exposure <- beta.outcome <- se.exposure <- se.outcome <- pval.exposure <- pval.outcome <- swap_unique_var_a <- NULL
  Ins <- data <- Ins
  Ins <- TwoSampleMR::format_data(Ins, type = "exposure", header = TRUE,
                     phenotype_col = "Phenotype", snp_col = "SNP", beta_col = "beta",
                     se_col = "se", eaf_col = "eaf", effect_allele_col = "effect_allele",
                     other_allele_col = "other_allele", pval_col = "pval")
# ao <- TwoSampleMR::available_outcomes(access_token=NULL)
  ids <- Ids
  outcome_dat <- TwoSampleMR::extract_outcome_data(snps = with(Ins,SNP), outcomes = ids)
  harmonise <- TwoSampleMR::harmonise_data(exposure_dat = Ins, outcome_dat = outcome_dat)
  if (reverse) harmonise <- subset(within(harmonise,
  {
    swap(exposure,outcome)
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
  result <- heterogeneity <- pleiotropy <- single <- NULL
  try(result <- TwoSampleMR::mr(harmonise, method_list=c("mr_wald_ratio", "mr_ivw"))) # main MR analysis
  heterogeneity <- TwoSampleMR::mr_heterogeneity(harmonise) # heterogeneity test across instruments
  pleiotropy <- TwoSampleMR::mr_pleiotropy_test(harmonise) # MR-Egger intercept test
  try(single <- TwoSampleMR::mr_singlesnp(harmonise)) #single SNP MR using Wald ratio
# result <- within(result,outcome <- sub(" [|]* id:ieu-a-[0-9]*| [|]* id:ukb-a-[0-9]*", "\\1", outcome, perl = TRUE))
  ext <- ".txt"
  invisible(lapply(c("harmonise","result","heterogeneity","pleiotropy","single"), function(x) {
                   v <- lapply(x, function(x) tryCatch(get(x), error=function(e) NULL))[[1]]
                   if (!is.null(v)) write.table(format(v,digits=3),file=paste0(prefix,"-",x,ext),quote=FALSE,row.names=FALSE,sep="\t")
            })
  )
}

uniprot2ids <- function(uniprotid="ACC+ID",to,query)
{
  rt <- find.package("pQTLtools")
  f <- file.path(rt ,"python","uniprot2ids.py")
  reticulate::source_python(f)
  invisible(uniprot2ids(uniprotid,to,query))
}

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

run_coloc <- function(eqtl_sumstats, gwas_sumstats)
{
  eQTL_dataset <- list(beta = eqtl_sumstats$beta,
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
  coloc_res <- coloc::coloc.abf(dataset1 = eQTL_dataset, dataset2 = gwas_dataset,p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
  res_formatted <- dplyr::as_tibble(t(as.data.frame(coloc_res$summary)))
  return(res_formatted)
}

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
