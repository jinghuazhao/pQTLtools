INF <- Sys.getenv("INF")
tsv <- paste(INF,"docs","hg19Tables.tsv",sep="/")
hg19Tables <- subset(read.delim(tsv,as.is=TRUE),!grepl("hap",X.chrom)&!grepl("Un",X.chrom)&!grepl("random",X.chrom)&!grepl(";",geneName)&geneName!="")
vars <- setdiff(names(hg19Tables),
                c("protFullNames","protShortNames","protAltFullNames","protAltShortNames","refSeq","refSeqProt","entrezGene","ensTrans","ensProt"))
hg19Tables <- hg19Tables[vars]
save(hg19Tables,file="hg19Tables.rda",compress='xz')

