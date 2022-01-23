#' A protein Quantitative Trait Locus toolkit
#'
#' The seeds collection of data and utilties for (pQTL) analysis. At
#' this early stage, this repository collects information on a number of
#' protein panels, linking function for cis/trans classification, 2D
#' manhattan plots, 3D-plotly plots, forest plots among others availale
#' from R/gap; query results on genes, regions, and SNPs via
#' PhenoScanner, adding functionality to check for replication across
#' platforms and aspects of protein-related analysis such as
#' pQTL-Mendelian Randomization via TwoSampleMR, linkage through UniProt
#' IDs to other resources.
#'
#' @details
#' \tabular{ll}{
#' Package: \tab pQTLtools\cr
#' Title: \tab pQTL Tools\cr
#' Author: \tab Jing Hua Zhao\cr
#' Maintainer: \tab Jing hua Zhao <jinghuazhao@hotmail.com>\cr
#' License: MIT + file LICENSE
#' URL: \tab https://github.com/jinghuazhao/pQTLtools\cr
#' BugReports: \tab https://github.com/jinghuazhao/pQTLtools/issues\cr
#' }
#'
#' @md
#'
#' ## A summary of datasets and functions
#'
#' Objects             |    Description
#' --------------------|-----------------------------------------
#' **Datasets**        |    
#' biomaRt             |    Curated data from biomaRt
#' caprion             |    Caprion panel
#' hg19                |    Curated data from Bioconductor
#' hg19Tables          |    Curated data from UCSC genome browser
#' inf1                |    Olink/INF panel
#' Olink_NGS           |    Olink/NGS panels
#' Olink_qPCR          |    Olink/qPCR panels
#' SomaLogic160410     |    SomaLogic panel
#' SomaScanV4.1        |    SomaScan v4.1 panel
#' st4                 |    ST4 of the INTERVAL SomaLogic paper
#' st6                 |    ST6 of the INTERVAL SomaLogic paper
#' st18                |    ST18 of the INTERVAL SomaLogic paper
#' swath_ms            |    SWATH-MS panel
#' **eQTL/GWAS**       |
#' get.prop.below.LLOD  |   Limit of detection analysis
#' import_eQTLCatalogue |   Import eQTL Catalogue
#' import_OpenGWAS      |   Import OpenGWAS
#' make_ExpressionSet   |   A call to ExpressionSet class
#' run_coloc            |   Colocalisation analysis
#' **MR analysis**      |
#' pqtlMR               |   Bidirectional pQTL-MR analysis
#' run_TwoSampleMR      |   A generic wrapper for TwoSampleMR analysis
#' **PhenoScanner Utilities** |
#' genequeries          |   phenoscanner genequeries in batches
#' regionqueries        |   phenoscanner regionqueries in batches
#' snpqueries           |   phenoscanner snpqueries in batches
#' **UniProt API**      |
#' uniprot2ids          |   UniProt ID to others
#'
#' Some generic description for the datasets are as follows.
#'
#' \describe{
#'   \item{chr}{chromosome}
#'   \item{start}{Start position}
#'   \item{end}{End position}
#'   \item{gene}{Gene name}
#'   \item{UniProt}{UniProt ID}
#' }
#' @examples
#' \dontrun{
#' # datasets
#' head(biomaRt)
#'
#' # Olink-SomaLogic panel overlap
#' p <- list(setdiff(inf1$uniprot,"P23560"),
#'           setdiff(SomaLogic160410$UniProt[!is.na(SomaLogic160410$UniProt)],"P23560"))
#' cnames <- c("INF1","SomaLogic")
#' VennDiagram::venn.diagram(x = p, category.names=cnames,
#'                           filename='os.png', imagetype="png", output=TRUE)
#' m <- merge(inf1,SomaLogic160410,by.x="uniprot",by.y="UniProt")
#' u <- setdiff(with(m,unique(uniprot)),"P23560")
#' options(width=220)
#' o <- subset(inf1,uniprot %in% u)
#' dim(o)
#' o
#' vars <- c("UniProt","chr","start","end","extGene","Target","TargetFullName")
#' s <- subset(SomaLogic160410[vars], UniProt %in% u)
#' dim(s)
#' us <- s[!duplicated(s),]
#' dim(us)
#' us
#'
#' # SCALLOP/INF1
#' INF <- Sys.getenv("INF")
#' INF1_merge <- merge(inf1,
#'                     read.delim(file.path(INF,"work","INF1.merge-rsid"),as.is=TRUE),
#'                     by="prot")
#' INF1_uniprot <- unique(with(INF1_merge,uniprot))
#'
#' # INTERVAL SomaLogic at box
#' HOME <- Sys.getenv("HOME")
#' box <- read.delim(file.path(HOME,"SomaLogic","doc","INTERVAL-box.tsv"),as.is=TRUE)
#' box_INF1 <- subset(box,UniProt %in% INF1_uniprot)
#' box_uniprot <- setdiff(unique(with(box_INF1,UniProt)),"P23560")
#' setdiff(INF1_uniprot,box_uniprot)
#'
#' # Phenoscanner database
#' ps <- merge(subset(read.delim(file.path(INF,"work","pQTL_2018.txt.gz"),as.is=TRUE),
#'             pmid==29875488),
#'             box,by.x="trait",by.y="TargetFullName")
#' z <- subset(ps,UniProtgwas %in% INF1_uniprot & p<=1.5e-11)
#'
#' # ST4 on Nature
#' st4regions <- subset(st4, UniProt %in% INF1_uniprot)
#' unique_uniprot_list <- setdiff(intersect(st4$UniProt,inf1$uniprot),"P23560")
#' subset(INF1_merge,uniprot %in% unique_uniprot_list)
#' }
#' @author Jing Hua Zhao in collaboration with other colleagues
#' @keywords internal

"_PACKAGE"
