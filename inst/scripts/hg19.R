library(OrganismDbi)
gd <- list(join1 = c(GO.db="GOID", org.Hs.eg.db="GO"), 
           join2 = c(org.Hs.eg.db = "ENTREZID", TxDb.Hsapiens.UCSC.hg19.knownGene = "GENEID"))
destination <- tempfile()
dir.create(destination)
makeOrganismPackage(pkgname = "Homo.sapiens.hg19", graphData = gd, 
                    organism = "Homo sapiens", version = "1.0.0", 
                    maintainer = "Maintainer<maintainer@email>", 
                    author = "Author Name", destDir = destination, 
                    license = "Artistic-2.0")
install.packages(paste(destination,"Homo.sapiens", sep="/"), repos = NULL, type="source")
library(Homo.sapiens.hg19)
Homo.sapiens.hg19
cols <- columns(Homo.sapiens.hg19)
cols
tx <- transcripts(Homo.sapiens.hg19, columns=c("UNIPROT","SYMBOL","ENSEMBL"))
chrs <- paste0("chr",c(1:22,'X','Y'))
hg19 <- subset(as.data.frame(tx),seqnames %in% chrs & !is.na(UNIPROT))
names(hg19)[1] <- "chr"

note <- function()
# Full genome sequence
{
  ip <- installed.packages()
  if (!"BSgenome.Hsapiens.UCSC.hg19" %in% rownames(ip)) BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
  library(BSgenome.Hsapiens.UCSC.hg19)
  Hsapiens
  library(BSgenome)
  gs <- getSeq(Hsapiens,seqnames(Hsapiens)[1:26])
}

save(hg19,file='hg19.rda',compress='xz')
