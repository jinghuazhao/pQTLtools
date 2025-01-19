#' QTL lookup
#'
#' This function takes MR results (involving pQTL/QTL) to look up QTLs in trait GWASs given a P value and
#' linkage disequilibrium (LD) cutoffs.
#' The rationale is that a pQTL may not necessarily be in strong LD with QTL but some other independent signal
#' in the same region. The interest then lands on association signals below a given P value (p_threshold) and above
#' a given LD (r2_threshold) thresholds.
#'
#' @param d directory where `dat` (below) is held.
#' @param dat MR results (`protein`, `id`, `pqtl`,`p`, `qtl`, `p_qtl`) whose `proxy`, `p_proxy` and `rsq` variables will be updated.
#' @param panel reference panel.
#' @param p_threshold cutoff of QTL association from trait GWASs.
#' @param r2_threshold cutoff of r^2.
#' @param pop reference population from 1000Genomes for LD calculation.
#' @param plink_bin PLINK executable file whose binary files are indicated in `bfile` variable of `dat`.
#' @param r when specified, the LD(r) is output.
#' @param r2 when specified, the LD(r^2) is output.
#' @param xlsx a non-null specification indicates name of an output Excel workbook.
#' 
#' @return A data.frame containing the looked up loci.
#' @export
#' @examples
#' \dontrun{
#' options(width=200)
#' INF <- Sys.getenv("INF")
#' suppressMessages(library(dplyr))
#' d <- file.path(INF,"mr","gsmr","trait")
#' inf1 <- select(gap.datasets::inf1,prot,target.short)
#' gsmr_efo <- read.delim(file.path(INF,"mr","gsmr","gsmr-efo.txt")) %>%
#'             left_join(inf1,by=c('protein'='target.short')) %>%
#'             mutate(file_gwas=paste(prot,id,"rsid.txt",sep="-"),
#'                    bfile=file.path(INF,"INTERVAL","per_chr",
#'                                    paste0("interval.imputed.olink.chr_",chr)),
#'                    proxy=NA,p_proxy=NA,rsq=NA)
#' proxies <- qtl_lookup(d,gsmr_efo,plink_bin="/rds/user/jhz22/hpc-work/bin/plink",
#'                       xlsx=file.path(INF,"mr","gsmr","r2_INTERVAL.xlsx")) %>%
#'            select(protein,id,Disease,fdr,pqtl,p,qtl,p_qtl,proxy,p_proxy,rsq)
#' write.table(proxies,file=file.path(INF,"mr","gsmr","r2_INTERVAL.tsv"),
#'             row.names=FALSE,quote=FALSE,sep="\t")
#' }

qtl_lookup <- function(d,dat,panel="1000Genomes",p_threshold=1e-3,r2_threshold=0.8,pop="EUR",
                       plink_bin=NULL,r=NULL,r2=NULL,xlsx=NULL)
{
  id <- protein <- p <- SNP <- prot <- Disease <- protein <- fdr <- qtl <- p_qtl <- rsq <- NULL
  for(i in 1:nrow(dat))
  {
     z <- dplyr::slice(dat,i)
     pqtl <- z[["pqtl"]]
     cat(z[["file_gwas"]],"\n")
     gwas <- read.table(file.path(d,basename(z[["file_gwas"]])),header=TRUE) %>%
             dplyr::arrange(p)
     h <-  dplyr::filter(gwas,p<=p_threshold)
     panel_snps <- c(pqtl,dplyr::pull(h,SNP))
     if (panel=="1000Genomes")
     {
        if (length(panel_snps)>500) stop("too many SNPs -- put a more stringent ptreshold")
        xx <- ieugwasr::ld_matrix(panel_snps,pop="EUR")
     } else
     xx <- ieugwasr::ld_matrix(panel_snps,with_alleles=TRUE,pop=pop,bfile=z[["bfile"]],plink_bin=plink_bin)
     cn <- colnames(xx)
     inside <- pqtl==gsub("_[A-Z]*","",cn)
     nn <- c(cn[inside],cn[!inside])
     if (length(r)==1) r_mat <- xx else r_mat <- xx[nn,nn]
     if(!is.null(r)) write.table(r_mat,file=paste(prot,id,Disease,pqtl,"r.txt",sep="-"))
     r2_mat <- r_mat^2
     if(!is.null(r2)) write.table(r2_mat,file=paste(prot,id,Disease,pqtl,"r2.txt",sep="-"))
     colnames(r2_mat) <- gsub("_[A-Z]*","",colnames(r2_mat))
     rownames(r2_mat) <- gsub("_[A-Z]*","",rownames(r2_mat))
     snps <- intersect(pull(h,SNP),colnames(r2_mat))
     success <- 1
     if (length(snps)<1) next else if (length(snps)==1) dat[i,c("proxy","p_proxy","rsq")] <- c(z[c("qtl","p_qtl")],1)
     else while(length(snps)>1)
     {
       success <- -1
       proxy <- snps[1]
       snps <- setdiff(snps,proxy)
       r2_i <- r2_mat[z[["pqtl"]],proxy]
       p_proxy <- filter(gwas,SNP==proxy) %>%
                  slice(which.min(p)) %>%
                  pull(p)
       cat("Same locus",i,z[["protein"]],z[["id"]],z[["Disease"]],z[["pqtl"]],z[["qtl"]],proxy,r2_i,z[["p_qtl"]],"\n",sep="\t")
       if(!is.null(r2_i)&!is.na(r2_i)) if(r2_i>r2_threshold) {success <-1; break}
     }
     if (success==-1)  while(length(snps)>1)
     {
       success <- -1
       proxy <- snps[1]
       snps <- setdiff(snps,proxy)
       r2_i <- r2_mat[z[["pqtl"]],proxy]
       p_proxy <- filter(gwas,SNP==proxy) %>%
                  slice(which.min(p)) %>%
                  pull(p)
       cat("Independent locus",i,z[["protein"]],z[["id"]],z[["Disease"]],z[["pqtl"]],z[["qtl"]],proxy,r2_i,z[["p_qtl"]],"\n",sep="\t")
       if(!is.null(r2_i)&!is.na(r2_i)) if(r2_i<r2_threshold) {success <-1; break}
     }
     if (success==1)
     {
       dat[i,"proxy"] <- proxy
       dat[i,"p_proxy"] <- p_proxy
       dat[i,"rsq"] <- r2_i
     }
  }
  if (!is.null(xlsx))
  {
    wb <- openxlsx::createWorkbook(xlsx)
    hs <- openxlsx::createStyle(textDecoration="BOLD", fontColour="#FFFFFF", fontSize=12, fontName="Arial Narrow", fgFill="#4F80BD")
    proxies <- dplyr::select(dat,protein,id,Disease,fdr,pqtl,p,qtl,p_qtl,proxy,p_proxy,rsq)
    for (sheet in "proxies")
    {
      openxlsx::addWorksheet(wb,sheet,zoom=150)
      openxlsx::writeData(wb,sheet,sheet,xy=c(1,1),headerStyle=openxlsx::createStyle(textDecoration="BOLD",
                          fontColour="#FFFFFF", fontSize=14, fontName="Arial Narrow", fgFill="#4F80BD"))
      body <- get(sheet)
      openxlsx::writeDataTable(wb, sheet, body, xy=c(1,2), headerStyle=hs, firstColumn=TRUE, tableStyle="TableStyleMedium2")
      openxlsx::freezePane(wb, sheet, firstActiveCol=3, firstActiveRow=3)
      width_vec <- apply(body, 2, function(x) max(nchar(as.character(x))+2, na.rm=TRUE))
    # width_vec_header <- nchar(colnames(body))+2
      openxlsx::setColWidths(wb, sheet, cols = 1:ncol(body), widths = width_vec)
      openxlsx::writeData(wb, sheet, tail(body,1), xy=c(1, nrow(body)+2), colNames=FALSE, borders="rows", borderStyle="thick")
    }
    openxlsx::saveWorkbook(wb, file=xlsx, overwrite=TRUE)
  }
  dat
}
