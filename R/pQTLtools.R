#' A reference manual
#'
#' @docType package
#' @name pQTLtools
#' @aliases pQTLtools-package
#'
#' @details
#' Available data and functions are listed in the following table.
#'
#' Objects                  |   Description
#' -------------------------|---------------------------------------------
#' **GWAS**                 |   &nbsp;
#' [`novelty_check`]        |   Locus novelty check
#' [`peptideAssociationPlot`]   | peptide association plot
#' [`peptideMapping`]       |   peptide-to-protein mapping
#' [`turboman`]             |   Manhattan plots
#' [`turboqq`]              |   QQ plots
#'  &nbsp;                  |   &nbsp;
#' **pGWAS**                |   &nbsp;
#' [`csq`]                  |   Variant consequence
#' [`protein_altering_variants`] | Protein Altering Variants (PAVs)
#' **Expression analysis**  |   &nbsp;
#' [`get.prop.below.LLOD`]  |   Limit of detection analysis
#' [`import_eQTLCatalogue`] |   Import eQTL Catalogue
#' [`make_ExpressionSet`]   |   A call to ExpressionSet class
#' [`run_coloc`]            |   Colocalisation analysis
#'  &nbsp;                  |   &nbsp;
#'  **MR analysis**         |   &nbsp;
#' [`import_OpenGWAS`]      |   Import OpenGWAS
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
#' [`gap::qtlFinder`]           |    Distance-based signal identification
#' [`gap::qtl2dplot`]           |    2D QTL plot
#' [`gap::qtl2dplotly`]         |    2D QTL plotly
#' [`gap::qtl3dplotly`]         |    3D QTL plotly
#'
#' @section Usage:
#' Vignettes on package usage:
#' - An Overview of pQTLtools. `vignette("pQTLtools")`.
#' - Bioconductor Notes. `vignette("bioconductor")`.
#' - ExpressionSet usage. `vignette("es")`.
#' - LocusZoom.js. `vignette("LocusZoom.js")`.
#' - snakemake showcases. `vignette("snakemake")`.
#' - Spectrum data analysis. `vignette("spectrum")`.
#' - SCALLOP-INF scripts. `vignette("SCALLOP-INF")`.
#' @md
#'
#' @import dplyr ggplot2 pQTLdata
#' @importFrom graphics axis legend lines mtext par points polygon segments text title
#' @importFrom grDevices colors xy.coords
#' @importFrom stats complete.cases median na.omit pnorm qbeta qchisq qnorm setNames
#' @importFrom utils head read.table tail write.table
#' @importFrom Rdpack reprompt
#'
#' @author Jing Hua Zhao in collaboration with other colleagues
#' @keywords internal

"_PACKAGE"
