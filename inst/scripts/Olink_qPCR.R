Lines <- 36
INF <- Sys.getenv("INF")
r <- readLines(paste(INF,"doc","Olink.R",sep='/'),n=Lines)
write(r,file=paste(Lines))
source(paste(Lines))
unlink(paste(Lines))
unlink("inf1.csv")
unlink("inf2.csv")

panels <- list()
for(i in 1:length(tabs))
{
  panel <- gsub(" |-", "_", tabs[i])
  target_uniprot <- get(panel)[1:2]
  names(target_uniprot) <- c("Target","UniProt")
  panels[[i]] <- data.frame(Panel=panel,target_uniprot)
}
p12 <- do.call(rbind,panels)

library(biomaRt)

# hg19/GRCh37
# listMarts()
# listEnsemblArchives()
hg19 <- useMart(biomart= "ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host = "grch37.ensembl.org")
# listDatasets(hg19)
# listAttributes(hg19)
# listFilters(hg19)
hg19.bm <- getBM(attributes = c('uniprotswissprot', 'hgnc_symbol','chromosome_name', 'start_position', 'end_position'),
                 filters = 'uniprotswissprot',
                 values = unique(with(p12,UniProt)),
                 mart = hg19)
hg19.bed <- subset(hg19.bm,chromosome_name%in%c(1:22,'X','Y'))
names(hg19.bed) <- c("UniProt","gene","chr","start","end")
Olink_qPCR <- merge(p12,hg19.bed,by="UniProt",all=TRUE)
save(Olink_qPCR,file='Olink_qPCR.rda',compress='xz')

single_panel <- function()
  save(Cardiometabolic, 
       Cell_Regulation,
       CVD_II,
       CVD_III,
       Development,
       Immune_Response,
       Immuno_Oncology,
       Inflammation,
       Metabolism,
       Neurology,
       Oncology_II,
       Organ_Damage,
       file='Olink_qPCR_info.rda',compress='xz')
