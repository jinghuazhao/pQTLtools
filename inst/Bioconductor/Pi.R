# https://bioconductor.org/packages/release/bioc/html/Pi.html

library(Pi)
RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
data.file <- "http://galahad.well.ox.ac.uk/bigdata/Spondyloarthritis.txt"
data <- read.delim(data.file, header=TRUE, stringsAsFactors=FALSE)
knitr::kable(data[1:10,], digits=200, caption="", row.names=FALSE)

include.LD <- 'EUR'
LD.r2 <- 0.8
LD.customised <- NULL
significance.threshold <- 5e-8
distance.max <- 50000
decay.kernel <- "rapid"
decay.exponent <- 2
include.eQTL <- c("JKscience_TS2A","JKscience_TS2B","JKscience_TS3A","JKng_bcell","JKng_mono","JKnc_neutro", "GTEx_V4_Whole_Blood")
eQTL.customised <- NULL
include.HiC <- c("Monocytes","Neutrophils","Total_B_cells")
GR.SNP <- "dbSNP_GWAS"
GR.Gene <- "UCSC_knownGene"
cdf.function <- "empirical"
relative.importance <- c(1/3, 1/3, 1/3)
scoring.scheme <- 'max'
network <- "STRING_high"
network.customised <- NULL
weighted <- FALSE
normalise <- "laplacian"
normalise.affinity.matrix <- "none"
restart <- 0.75
parallel <- TRUE
multicores <- NULL
verbose <- TRUE

pNode <- xPierSNPs(data, include.LD=include.LD, LD.r2=LD.r2, significance.threshold=significance.threshold,
                   distance.max=distance.max, decay.kernel=decay.kernel, decay.exponent=decay.exponent,
                   include.eQTL=include.eQTL, include.HiC=include.HiC, GR.SNP=GR.SNP, GR.Gene=GR.Gene,
                   cdf.function=cdf.function, scoring.scheme=scoring.scheme, network=network, restart=restart, RData.location=RData.location)
write.table(pNode$priority, file="Genes_priority.txt", sep="\t", row.names=FALSE)
mp <- xPierManhattan(pNode, top=40, y.scale="sqrt", RData.location=RData.location)
print(mp)
png(file=file.path("figures","Pi.gene_manhattan.png"), height=480, width=480*2.2, res=96*1.3)
print(mp)
dev.off()
library(png)
library(grid)
library(XGR)
eTerm <- xPierPathways(pNode, priority.top=100, ontology="MsigdbC2CPall", RData.location=RData.location)
eTerm_nonred <- xEnrichConciser(eTerm)

# view the top pathways/terms
xEnrichViewer(eTerm_nonred)
Pathways_nonred <- xEnrichViewer(eTerm_nonred, top_num=length(eTerm_nonred$adjp), sortBy="fdr", details=TRUE)
output <- data.frame(term=rownames(Pathways_nonred), Pathways_nonred)
write.table(output, file="Pathways_priority.nonredundant.txt", sep="\t", row.names=FALSE)
bp <- xEnrichBarplot(eTerm_nonred, top_num="auto", displayBy="fdr", FDR.cutoff=1e-3, wrap.width=50, signature=FALSE)
bp
png(file=file.path("figures","Pi.pathway_barplot.png"), height=360, width=360*3, res=96)
print(bp)
dev.off()

# find maximum-scoring gene network with the desired node number=50
g <- xPierSubnet(pNode, priority.quantile=0.1, subnet.size=50, RData.location=RData.location)

pattern <- as.numeric(V(g)$priority)
zmax <- ceiling(quantile(pattern,0.75)*1000)/1000
xVisNet(g, pattern=pattern, vertex.shape="sphere", colormap="yr", newpage=FALSE, edge.arrow.size=0.3,
        vertex.label.color="blue", vertex.label.dist=0.35, vertex.label.font=2, zlim=c(0,zmax), signature=FALSE)

png(file=file.path("figures","Pi.network_vis.png"), height=480*2, width=480*2, res=96*1.5)
xVisNet(g, pattern=pattern, vertex.shape="sphere", colormap="yr", newpage=FALSE, vertex.label.color="blue",
        vertex.label.dist=0.35, vertex.label.font=2, zlim=c(0,zmax), signature=FALSE)
dev.off()
