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
#' outcome <- TwoSampleMR::extract_outcome_data(snps=exposure$SNP,outcomes="ebi-a-GCST007432")
#' harmonise <- TwoSampleMR::harmonise_data(clump,outcome)
#' prefix <- paste(prot,type,sep="-")
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
