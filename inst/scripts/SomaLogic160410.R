HOME <- Sys.getenv("HOME")

# from Cardio
SomaLogic <- read.delim(paste(HOME,"SomaLogic","doc","SOMALOGIC_Master_Table_160410_1129info.tsv",sep="/"),as.is=TRUE)
vars <- c("SOMAS_ID_round2", "UniProt", "Target", "TargetFullName",
          "chromosome_number","start_position","end_position","EntrezGeneSymbol","ensembl_gene_id","external_gene_name")
chrs <- c(paste(1:22),"X","Y")
library(reshape)
SomaLogic160410 <- rename(subset(SomaLogic[vars],chromosome_number%in%chrs),
                          c(chromosome_number="chr",start_position="start",end_position="end",
                            EntrezGeneSymbol="entGene",ensembl_gene_id="ensGene",external_gene_name="extGene",
                            SOMAS_ID_round2="SOMAMER_ID"))
out <- c("chr","start","end","ensGene","UniProt","entGene","TargetFullName","extGene","Target","UniProt","SOMAMER_ID")
save(SomaLogic160410, file="SomaLogic160410.rda", compress='xz')

# from Box
gwas <- read.csv(paste(HOME,"SomaLogic","doc","SOMALOGIC_GWAS_protein_info.csv",sep="/"),as.is=TRUE)
gs <- merge(gwas[c("SOMAMER_ID","Target","TargetFullName")],SomaLogic160410[setdiff(out,c("Target","TargetFullName"))],by="SOMAMER_ID")
ord <- with(gs,order(chr,start,end))
INTERVAL_gwas <- gs[ord,out]

# INTERVAL-box.tsv is available from the SomaLogic GitHub repository
write.table(INTERVAL_gwas,file="INTERVAL-box.tsv",quote=FALSE,row.names=FALSE,sep="\t")
