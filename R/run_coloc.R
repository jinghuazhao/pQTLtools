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
