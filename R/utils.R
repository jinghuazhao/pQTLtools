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

#' Colocalisation analysis
#'
#' The function takes eQTL and GWAS summary statistics for a colocalisation analysis..
#'
#' @md
#' @param eqtl_sumstats eQTL summary data.
#' @param gwas_sumstats GWAS summary data.
#' @export
#' @return Summary from `coloc.abf`.

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


blank_plot <- function(message)
{
   requireNamespace("ggplot2", quietly=TRUE)
   a <- b <- NA
   ggplot2::ggplot(data.frame(a=0,b=0,n=message)) +
   ggplot2::geom_text(ggplot2::aes(x=a,y=b,label=n)) +
   ggplot2::labs(x=NULL,y=NULL) +
   ggplot2::theme(axis.text=ggplot2::element_blank(), axis.ticks=ggplot2::element_blank())
}

mr_scatter_plot2 <- function (mr_results, dat, alpha=0.05)
{
   requireNamespace("ggplot2", quietly = TRUE)
   requireNamespace("plyr", quietly = TRUE)
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

mr_forest_plot2 <- function (singlesnp_results, exponentiate = FALSE, alpha = 0.05)
{
  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("plyr", quietly = TRUE)
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
          cowplot::theme_cowplot(12) +
          ggplot2::geom_vline(xintercept = xint, linetype = "dotted") +
          ggplot2::geom_errorbarh(ggplot2::aes(xmin = lo, xmax = up, size = as.factor(tot), colour = as.factor(tot)), height = 0) +
          ggplot2::geom_point(ggplot2::aes(colour = as.factor(tot))) +
          ggplot2::geom_hline(ggplot2::aes(yintercept = which(levels(SNP) %in% ""))) +
          ggplot2::scale_colour_manual(values = c("black", "red")) +
          ggplot2::scale_size_manual(values = c(0.3, 1)) +
          ggplot2::theme(legend.position = "none", axis.text.y = ggplot2::element_text(size = 8),
                         axis.ticks.y = ggplot2::element_line(size = 0), axis.title.x = ggplot2::element_text(size = 8)) +
          ggplot2::labs(y = "", x = paste0("MR effect size for\n'", with(d, exposure)[1], "' on '", with(d, outcome)[1], "'"))
      })
  res
}

mr_funnel_plot2 <- function (singlesnp_results)
{
  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("plyr", quietly = TRUE)
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

mr_leaveoneout_plot2 <- function (leaveoneout_results, alpha = 0.05)
{
  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("plyr", quietly = TRUE)
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
      cowplot::theme_cowplot(12) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dotted") +
      ggplot2::geom_errorbarh(ggplot2::aes(xmin = lo, xmax = up, size = as.factor(tot), colour = as.factor(tot)), height = 0) +
      ggplot2::geom_point(ggplot2::aes(colour = as.factor(tot))) +
      ggplot2::geom_hline(ggplot2::aes(yintercept = which(levels(SNP) %in% ""))) +
      ggplot2::scale_colour_manual(values = c("black", "red")) +
      ggplot2::scale_size_manual(values = c(0.3, 1)) +
      ggplot2::theme(legend.position = "none", axis.text.y = ggplot2::element_text(size = 8),
                     axis.ticks.y = ggplot2::element_line(size = 0), axis.title.x = ggplot2::element_text(size = 8)) +
      ggplot2::labs(y = "", x = paste0("MR leave-one-out sensitivity analysis for\n'", with(d, exposure)[1], "' on '", with(d, outcome)[1], "'"))
  })
  res
}

pqtlMR <- function(ivs, ids, mr_plot=FALSE, prefix="pQTL-combined-", reverse=FALSE)
{
  outcome_data <- extract_outcome_data(snps=with(ivs,SNP),outcomes=ids)
  harmonise <- harmonise_data(ivs,outcome_data)
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
# main MR analysis
  try(result <- TwoSampleMR::mr(harmonise,method_list=c("mr_wald_ratio", "mr_ivw")))
# single SNP MR using Wald ratio
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

run_TwoSampleMR <- function(TwoSampleMRinput, mr_plot="None", prefix="")
{
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
#' library(Biobase)
#' data(sample.ExpressionSet)
#' exampleSet <- sample.ExpressionSet
#' fData(exampleSet)
#' fData(exampleSet)$lod.max <- apply(exprs(exampleSet),1,quantile,runif(nrow(exampleSet)))
#' lod <- get.prop.below.LLOD(exampleSet)
#' x <- dplyr::arrange(fData(lod),desc(pc.belowLOD.new))
#' knitr::kable(head(lod))
#' plot(x[,2], main="Random quantile cut off", ylab="<lod%")

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

  requireNamespace("stringr")

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
