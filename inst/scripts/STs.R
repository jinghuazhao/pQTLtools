# INTERVAL SomaLogic results
dir <- 'https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-018-0175-2/MediaObjects/'
file <- '41586_2018_175_MOESM4_ESM.xlsx'
xlsx <- paste0(dir,file)
st4.1 <- openxlsx::read.xlsx(xlsx, sheet=4, colNames=TRUE, skipEmptyRows=TRUE,
                             cols=c(1:16,26:28), rows=c(5:1986))
st4.2 <- openxlsx::read.xlsx(xlsx, sheet=4, colNames=TRUE, skipEmptyRows=TRUE,
                             cols=c(17:25,29:31), rows=c(6:1986))
st4 <- cbind(st4.1,st4.2)
save(st4,file='st4.rda',compress='xz')
st5 <- openxlsx::read.xlsx(xlsx, sheet=5, colNames=TRUE, skipEmptyRows=TRUE,
                           cols=c(1:19), rows=c(3:2746))
st6.1 <- openxlsx::read.xlsx(xlsx, sheet=6, colNames=TRUE, skipEmptyRows=TRUE,
                           cols=c(1:11:13), rows=c(3:167))
st6.2 <- openxlsx::read.xlsx(xlsx, sheet=6, colNames=TRUE, skipEmptyRows=TRUE,
                           cols=c(14:20), rows=c(4:167))
st6 <- cbind(st6.1,st6.2)
save(st6,file='st6.rda',compress='xz')
replicates <- merge(st4[,c(1:12,26:28)],st6[,c(1:10,17:20)],
                    by=c("Locus.ID","UniProt","Chr","Pos","SOMAmer.ID"))
line_to_edit <- with(replicates,Locus.ID=="3_59"|Locus.ID=="5_29")
replicates[line_to_edit,"UniProt"] <- "P29460"

# However, it is unclear UniProts in ST6 were selected from which of the panels
INF <- Sys.getenv("INF")
INF1_merge <- merge(inf1,
                    read.delim(file.path(INF,"work","INF1.merge-rsid"),as.is=TRUE),
                    by="prot")
INF1_uniprot <- unique(with(INF1_merge,uniprot))
options(width=250)
st6_replicates <- subset(replicates,UniProt %in% INF1_uniprot)
table(subset(st6_replicates,UniProt %in% INF1_uniprot)$Replicates)

# side information on cvd2, cvd3, inf1
olink <- scan(file.path(INF,"doc","olink.prot.list.txt"),"")
olink_uniprot <- unlist(lapply(strsplit(olink,"___"),'[[',2))
dim(subset(replicates,UniProt %in% olink_uniprot))

z <- with(st6_replicates,{
  chr <- st6_replicates[["Chr"]]
  pos <- st6_replicates[["Pos"]]
  a1 <- st6_replicates[["Effect.Allele.(EA)"]]
  a2 <- st6_replicates[["Other.Allele.(OA)"]]
  cbind(UniProt,snpid=chr_pos_a1_a2(chr,pos,a1,a2))
})

doubles <- c("P29460","Q9NPF7","Q14213","Q8NEV9")
subset(inf1,uniprot %in% doubles)
write.table(merge(inf1,z,by.x="uniprot",by.y="UniProt")[c("snpid","prot","uniprot")],
            file="SomaLogic.id3",col.names=FALSE,row.names=FALSE,quote=FALSE)

st18 <- openxlsx::readWorkbook(xlsx,sheet="ST18 - proteins assayed",startRow=3)
save(st18,file='st18.rda',compress='xz')
st20 <- openxlsx::read.xlsx(xlsx, sheet=20, colNames=TRUE, skipEmptyRows=TRUE,
                            cols=c(1:10), rows=c(3:786))
ov <- intersect(unique(st20$UniProt),inf1$uniprot)
ovv <- subset(st20,UniProt%in%ov)
dim(ovv)
