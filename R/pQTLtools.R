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
#' Objects              \tab   Description\cr
#' \cr
#' \strong{eQTL/GWAS}\cr
#' \cr
#' `get.prop.below.LLOD`  \tab   Limit of detection analysis\cr
#' `import_eQTLCatalogue` \tab   Import eQTL Catalogue\cr
#' `import_OpenGWAS`      \tab   Import OpenGWAS\cr
#' `make_ExpressionSet`   \tab   A call to ExpressionSet class\cr
#' `novelty_check`        \tab   Locus novelty check\cr
#' `run_coloc`            \tab   Colocalisation analysis\cr
#' \cr
#' \strong{MR analysis}\cr
#' \cr
#' `pqtlMR`               \tab   Bidirectional pQTL-MR analysis\cr
#' `qtl_lookup`           \tab   QTL lookup\cr
#' `run_TwoSampleMR`      \tab   A generic wrapper for TwoSampleMR analysis\cr
#' \cr
#' \strong{PhenoScanner Utilities}\cr
#' \cr
#' `genequeries`          \tab   phenoscanner genequeries in batches\cr
#' `regionqueries`        \tab   phenoscanner regionqueries in batches\cr
#' `snpqueries`           \tab   phenoscanner snpqueries in batches\cr
#' \cr
#' \strong{UniProt API}   \cr
#' \cr
#' `uniprot2ids`          \tab   UniProt ID to others\cr
#' \cr
#' \strong{Functions in gap}\cr
#' \cr
#' [`gap::METAL_forestplot`]    \tab    Forest plots from metal analysis\cr
#' [`gap::ci2ms`]               \tab    Effect size and standard error from confidence interval\cr
#' [`gap::cis.vs.trans.classification`] \tab a cis/trans classifier\cr
#' [`gap::circos.cis.vs.trans.plot`] \tab circos plot of cis/trans classification\cr
#' [`gap::circos.mhtplot`]      \tab    circos Manhattan plot with gene annotation\cr
#' [`gap::circos.mhtplot2`]     \tab    Another circos Manhattan plot\cr
#' [`gap::cs`]                  \tab    Credible set\cr
#' [`gap::get_b_se`]            \tab    Get b and se from AF, n, and z\cr
#' [`gap::get_pve_se`]          \tab    Get pve and its standard error from n, z\cr
#' [`gap::get_sdy`]             \tab    Get sd(y) from AF, n, b, se\cr
#' [`gap::mr`]                  \tab    Mendelian randomization analysis\cr
#' [`gap::invnormal`]           \tab    Inverse normal transformation\cr
#' [`gap::log10p`]              \tab    log10(p) for a standard normal deviate\cr
#' [`gap::log10pvalue`]         \tab    log10(p) for a P value including its scientific format\cr
#' [`gap::logp`]                \tab    log(p) for a normal deviate\cr
#' [`gap::mhtplot.trunc`]       \tab    Truncated Manhattan plot\cr
#' [`gap::miamiplot2`]          \tab    Miami plot\cr
#' [`gap::mr_forestplot`]       \tab    Mendelian Randomization forest plot\cr
#' [`gap::qtlClassifier`]       \tab    A QTL cis/trans classifier\cr
#' [`gap::qtl2dplot`]           \tab    2D QTL plot\cr
#' [`gap::qtl2dplotly`]         \tab    2D QTL plotly\cr
#' [`gap::qtl3dplotly`]         \tab    3D QTL plotly
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
