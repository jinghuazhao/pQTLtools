#' A summary of functions
#'
#' It seeds collection of data and utilities for pQTL analysis. At this
#' early stage, the repository contains 1. Articles linking functions
#' for cis/trans classification, 2d Manhattan plots, 2d/3d-plotly plots,
#' forest plots among others available from
#' gap, <https://github.com/jinghuazhao/R/tree/master/gap>;
#' 2. Query on genes, regions, and SNPs via PhenoScanner,
#' <http://www.phenoscanner.medschl.cam.ac.uk/>, adding functionality
#' to check for replication across platforms; 3. Downstream analysis
#' such as colocalization, pQTL-Mendelian Randomization via TwoSampleMR,
#' <https://github.com/MRCIEU/TwoSampleMR>, linkage through UniProt IDs
#' to other resources; 4. Bioconductor notes and a showcase of snakemake
#' workflow.
#'
#' @details
#' Available functions are listed in the following table.
#'
#' \tabular{ll}{
#' Objects              &   Description\cr
#' \cr
#' \strong{eQTL/GWAS}\cr
#' \cr
#' `get.prop.below.LLOD`  &   Limit of detection analysis\cr
#' `import_eQTLCatalogue` &   Import eQTL Catalogue\cr
#' `import_OpenGWAS`      &   Import OpenGWAS\cr
#' `make_ExpressionSet`   &   A call to ExpressionSet class\cr
#' `novelty_check`        &   Locus novelty check\cr
#' `run_coloc`            &   Colocalisation analysis\cr
#' \cr
#' \strong{MR analysis}\cr
#' \cr
#' `pqtlMR`               &   Bidirectional pQTL-MR analysis\cr
#' `qtl_lookup`           &   QTL lookup\cr
#' `run_TwoSampleMR`      &   A generic wrapper for TwoSampleMR analysis\cr
#' \cr
#' \strong{PhenoScanner Utilities}\cr
#' \cr
#' `genequeries`          &   phenoscanner genequeries in batches\cr
#' `regionqueries`        &   phenoscanner regionqueries in batches\cr
#' `snpqueries`           &   phenoscanner snpqueries in batches\cr
#' \cr
#' \strong{UniProt API}   \cr
#' \cr
#' `uniprot2ids`          &   UniProt ID to others\cr
#' \cr
#' \strong{Functions in gap}\cr
#' \cr
#' [`gap::METAL_forestplot`]    &    Forest plots from metal analysis\cr
#' [`gap::ci2ms`]               &    Effect size and standard error from confidence interval\cr
#' [`gap::cis.vs.trans.classification`] & a cis/trans classifier\cr
#' [`gap::circos.cis.vs.trans.plot`] & circos plot of cis/trans classification\cr
#' [`gap::circos.mhtplot`]      &    circos Manhattan plot with gene annotation\cr
#' [`gap::circos.mhtplot2`]     &    Another circos Manhattan plot\cr
#' [`gap::cs`]                  &    Credible set\cr
#' [`gap::get_b_se`]            &    Get b and se from AF, n, and z\cr
#' [`gap::get_pve_se`]          &    Get pve and its standard error from n, z\cr
#' [`gap::get_sdy`]             &    Get sd(y) from AF, n, b, se\cr
#' [`gap::mr`]                  &    Mendelian randomization analysis\cr
#' [`gap::invnormal`]           &    Inverse normal transformation\cr
#' [`gap::log10p`]              &    log10(p) for a standard normal deviate\cr
#' [`gap::log10pvalue`]         &    log10(p) for a P value including its scientific format\cr
#' [`gap::logp`]                &    log(p) for a normal deviate\cr
#' [`gap::mhtplot.trunc`]       &    Truncated Manhattan plot\cr
#' [`gap::miamiplot2`]          &    Miami plot\cr
#' [`gap::mr_forestplot`]       &    Mendelian Randomization forest plot\cr
#' [`gap::qtlClassifier`]       &    A QTL cis/trans classifier\cr
#' [`gap::qtl2dplot`]           &    2D QTL plot\cr
#' [`gap::qtl2dplotly`]         &    2D QTL plotly\cr
#' [`gap::qtl3dplotly`]         &    3D QTL plotly
#' }
#'
#' @section Usage:
#' Vignettes on package usage:
#' - An Overview of pQTLtools. `vignette("pQTLtools")`.
#' - Bioconductor Notes. `vignette("bioconductor")`.
#' - ExpressionSet usage. `vignette("es")`.
#' - snakemake showcases. `vignette("snakemake")`.
#' - SCALLOP-INF scripts. `vignette("SCALLOP-INF")`.
#' @md
#'
#' @docType package
#' @name pQTLtools
#' @aliases pQTLtools-package
#'
#' @import dplyr ggplot2 pQTLdata
#' @importFrom stats na.omit qnorm
#' @importFrom utils write.table
#' @importFrom utils read.table tail
#' @importFrom Rdpack reprompt
#'
#' @author Jing Hua Zhao in collaboration with other colleagues
#' @keywords internal

"_PACKAGE"
