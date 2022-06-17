#' A summary of datasets and functions
#'
#' \tabular{ll}{
#' Objects             \tab    Description\cr
#' \strong{Datasets}   \tab    \cr
#' biomaRt             \tab    Curated data from biomaRt\cr
#' caprion             \tab    Caprion panel\cr
#' hg19                \tab    Curated data from Bioconductor\cr
#' hg19Tables          \tab    Curated data from UCSC genome browser\cr
#' inf1                \tab    Olink/INF panel\cr
#' Olink_NGS           \tab    Olink/NGS panels\cr
#' Olink_qPCR          \tab    Olink/qPCR panels\cr
#' SomaLogic160410     \tab    SomaLogic panel\cr
#' SomaScanV4.1        \tab    SomaScan v4.1 panel\cr
#' st4                 \tab    ST4 of the INTERVAL SomaLogic paper\cr
#' st6                 \tab    ST6 of the INTERVAL SomaLogic paper\cr
#' st18                \tab    ST18 of the INTERVAL SomaLogic paper\cr
#' swath_ms            \tab    SWATH-MS panel\cr
#' \strong{eQTL/GWAS}  \tab    \cr
#' get.prop.below.LLOD  \tab   Limit of detection analysis\cr
#' import_eQTLCatalogue \tab   Import eQTL Catalogue\cr
#' import_OpenGWAS      \tab   Import OpenGWAS\cr
#' make_ExpressionSet   \tab   A call to ExpressionSet class\cr
#' run_coloc            \tab   Colocalisation analysis\cr
#' \strong{MR analysis} \tab   \cr
#' pqtlMR               \tab   Bidirectional pQTL-MR analysis\cr
#' run_TwoSampleMR      \tab   A generic wrapper for TwoSampleMR analysis\cr
#' \strong{PhenoScanner Utilities} \tab \cr
#' genequeries          \tab   phenoscanner genequeries in batches\cr
#' regionqueries        \tab   phenoscanner regionqueries in batches\cr
#' snpqueries           \tab   phenoscanner snpqueries in batches\cr
#' \strong{UniProt API} \tab   \cr
#' uniprot2ids          \tab   UniProt ID to others\cr
#' }
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
#'
#' @import ggplot2
#' @importFrom Biobase annotatedDataFrameFrom MIAME ExpressionSet exprs fData
#' @importFrom GenomicRanges GRanges width tile
#' @importFrom IRanges IRanges
#' @importFrom VariantAnnotation VcfFile vcfFields ScanVcfParam readVcf
#' @importFrom coloc coloc.abf
#' @importFrom cowplot theme_cowplot
#' @importFrom dplyr as_tibble filter select distinct mutate group_by ungroup
#' @importFrom gwasvcf set_bcftools query_gwas vcf_to_granges
#' @importFrom ieugwasr check_access_token
#' @importFrom phenoscanner phenoscanner
#' @importFrom plyr ldply
#' @importFrom reticulate source_python
#' @importFrom rtracklayer import.chain liftOver
#' @importFrom seqminer tabix.read.table
#' @importFrom TwoSampleMR mr_egger_regression mr_egger_regression_bootstrap mr mr_singlesnp mr_scatter_plot mr_forest_plot mr_funnel_plot mr_heterogeneity mr_pleiotropy_test mr_leaveoneout mr_leaveoneout_plot
#'
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
