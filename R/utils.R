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

extract_outcome_data.args <- function(snps=NULL, outcomes=NULL, proxies=TRUE, rsq=0.8, align_alleles=1,
                                      palindromes=1, maf_threshold=0.3, access_token=ieugwasr::check_access_token(),
                                      splitsize=10000, proxy_splitsize=500)
  invisible(list(snps=snps,outcomes=outcomes,proxies=proxies,rsq=rsq,align_alleles=align_alleles,
                           palindromes=palindromes,maf_threshold=maf_threshold,access_token=access_token,
                           splitsize=splitsize,proxy_splitsize=proxy_splitsize))

clump_data.args <- function(dat=NULL, clump_kb = 10000, clump_r2 = 0.001, clump_p1 = 1, clump_p2 = 1, pop = "EUR")
  invisible(list(dat=dat,clump_kb=clump_kb,clump_r2=clump_r2,clump_p1=clump_p1,clump_p2=clump_p2,pop=pop))

harmonise_data.args <- function(exposure_dat=NULL, outcome_dat=NULL, action = 2)
  invisible(list(exporesure_dat=exposure_dat,outcome_dat=outcome_dat,action=action))

swap <- function(x,y)
   eval(parse(text = paste("swap_unique_var_a <-", substitute(x), ";",
   substitute(x), "<-", substitute(y), ";",
   substitute(y), "<-swap_unique_var_a")), envir=parent.frame())

pqtlMR <- function(pqtlMRinput, plot=TRUE, prefix="pQTL-combined-",reverse=FALSE)
{
  Ins <- with(pqtlMRinput,Ins)
  Ids <- with(pqtlMRinput,Ids)
  harmonise <- with(pqtlMRinput,harmonise)
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
  if (is.null(Ins) | is.null(Ids) | is.null(harmonise))
     stop("Missing element(s) in the input")
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
# main MR analysis
  try(result <- TwoSampleMR::mr(harmonise,method_list=c("mr_wald_ratio", "mr_ivw")))
# single SNP MR using Wald ratio
  try(single <- TwoSampleMR::mr_singlesnp(harmonise))
  invisible(lapply(c("harmonise","result","single"), function(x) if (exists(x))
            write.table(format(get(x),digits=3),file=paste0(prefix,x,".txt"),
                        quote=FALSE,row.names=FALSE,sep="\t")))
  if (plot)
  {
    scatter <- TwoSampleMR::mr_scatter_plot(result, harmonise)
    forest <- TwoSampleMR::mr_forest_plot(single)
    funnel <- TwoSampleMR::mr_funnel_plot(single)
    invisible(sapply(c(scatter,forest,funnel), function(x) print(x)))
  }
}

run_TwoSampleMR <- function(TwoSampleMRinput, plot=TRUE, prefix="")
{
  exposure <- with(TwoSampleMRinput, exposure)
  outcome <- with(TwoSampleMRinput, outcome)
  clump <- with(TwoSampleMRinput, clump)
  harmonise <- with(TwoSampleMRinput, harmonise)
  if (is.null(exposure) | is.null(outcome) | is.null(clump) | is.null(harmonise))
     stop("Missing element(s) in the input")
  TwoSampleMR::directionality_test(harmonise)
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
  if (plot)
  {
    scatter <- TwoSampleMR::mr_scatter_plot(result, harmonise)
    forest <- TwoSampleMR::mr_forest_plot(single)
    funnel <- TwoSampleMR::mr_funnel_plot(single)
    leaveoneout <- TwoSampleMR::mr_leaveoneout_plot(loo)
    invisible(sapply(c(scatter,forest,funnel,leaveoneout), function(x) print(x)))
  }
}
