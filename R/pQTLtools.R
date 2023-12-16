#' A summary of functions
#'
#' It seeds collection of data and utilities for pQTL analysis. At this
#' early stage, the repository contains 1. Articles linking functions
#' for cis/trans classification, pQTL-gene plot, 2d/3d-plotly plots,
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
#' Objects                  |   Description
#' -------------------------|---------------------------------------------
#' **eQTL/GWAS**            |   &nbsp;
#' [`get.prop.below.LLOD`]  |   Limit of detection analysis
#' [`import_eQTLCatalogue`] |   Import eQTL Catalogue
#' [`import_OpenGWAS`]      |   Import OpenGWAS
#' [`make_ExpressionSet`]   |   A call to ExpressionSet class
#' [`novelty_check`]        |   Locus novelty check
#' [`run_coloc`]            |   Colocalisation analysis
#' [`turboman`]             |   Manhattan plots for impatient people
#' [`turboqq`]              |   QQ plots for impatient people
#'  &nbsp;                  |   &nbsp;
#'  **MR analysis**         |   &nbsp;
#' [`pqtlMR`]               |   Bidirectional pQTL-MR analysis
#' [`qtl_lookup`]           |   QTL lookup
#' [`run_TwoSampleMR`]      |   A generic wrapper for TwoSampleMR analysis
#'  &nbsp;                  |   &nbsp;
#'  **PhenoScanner Utilities** | &nbsp;
#' [`genequeries`]          |   phenoscanner genequeries in batches
#' [`regionqueries`]        |   phenoscanner regionqueries in batches
#' [`snpqueries`]           |   phenoscanner snpqueries in batches
#'  &nbsp;                  |   &nbsp;
#' **UniProt API**          |   &nbsp;
#' [`uniprot2ids`]          |   UniProt ID to others
#' &nbsp;                   |   &nbsp;
#' **Functions in gap**     |   &nbsp;
#' [`gap::METAL_forestplot`]    |    Forest plots from metal analysis
#' [`gap::ci2ms`]               |    Effect size and standard error from confidence interval
#' [`gap::cis.vs.trans.classification`] | a cis/trans classifier
#' [`gap::circos.cis.vs.trans.plot`] | circos plot of cis/trans classification
#' [`gap::circos.mhtplot`]      |    circos Manhattan plot with gene annotation
#' [`gap::circos.mhtplot2`]     |    Another circos Manhattan plot
#' [`gap::cs`]                  |    Credible set
#' [`gap::get_b_se`]            |    Get b and se from AF, n, and z
#' [`gap::get_pve_se`]          |    Get pve and its standard error from n, z
#' [`gap::get_sdy`]             |    Get sd(y) from AF, n, b, se
#' [`gap::mr`]                  |    Mendelian randomization analysis
#' [`gap::invnormal`]           |    Inverse normal transformation
#' [`gap::log10p`]              |    log10(p) for a standard normal deviate
#' [`gap::log10pvalue`]         |    log10(p) for a P value including its scientific format
#' [`gap::logp`]                |    log(p) for a normal deviate
#' [`gap::mhtplot.trunc`]       |    Truncated Manhattan plot
#' [`gap::miamiplot2`]          |    Miami plot
#' [`gap::mr_forestplot`]       |    Mendelian Randomization forest plot
#' [`gap::pvalue`]              |    P value for a normal deviate
#' [`gap::qtlClassifier`]       |    A QTL cis/trans classifier
#' [`gap::qtl2dplot`]           |    2D QTL plot
#' [`gap::qtl2dplotly`]         |    2D QTL plotly
#' [`gap::qtl3dplotly`]         |    3D QTL plotly
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
