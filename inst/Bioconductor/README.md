## Transcript databases

```r
options(width=200)

suppressMessages(library(AnnotationDbi))
suppressMessages(library(org.Hs.eg.db))
columns(org.Hs.eg.db)
keyref <- keys(org.Hs.eg.db, keytype="ENTREZID")
symbol_uniprot <- select(org.Hs.eg.db,keys=keyref,columns = c("SYMBOL","UNIPROT"))
subset(symbol_uniprot,SYMBOL=="MC4R")

suppressMessages(library(EnsDb.Hsapiens.v86))
x <- EnsDb.Hsapiens.v86
listColumns(x, "protein", skip.keys=TRUE)
listGenebiotypes(x)
listTxbiotypes(x)
listTables(x)
metadata(x)
organism(x)
returnFilterColumns(x)
seqinfo(x)
seqlevels(x)
updateEnsDb(x)

suppressMessages(library(ensembldb))
genes(x, columns=c("gene_name"), filter=list(SeqNameFilter("X"), GeneBiotypeFilter("protein_coding")))
transcripts(x, columns=listColumns(x, "tx"), filter = AnnotationFilterList(), order.type = "asc", return.type = "GRanges")

library(RMariaDB)
library(GenomicFeatures)
txdbEnsemblGRCh38 <- makeTxDbFromEnsembl(organism="Homo sapiens", release=98)
txdb <- as.list(txdbEnsemblGRCh38)
lapply(txdb,head)

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

library(INSPEcT)
liverExprs <- quantifyExpressionsFromBWs(txdb = txdb,BWfiles=,experimentalDesign=)
```

## Pathway and enrichment analysis

```r
library(graphite)
reactome <- pathways("hsapiens", "reactome")
kegg <- pathways("hsapiens","kegg")
pharmgkb <- pathways("hsapiens","pharmgkb")
nodes(kegg)
suppressMessages(library(plyr))
kegg_t2g <- ldply(lapply(kegg, nodes), data.frame)
names(kegg_t2g) <- c("gs_name", "gene_symbol")
suppress(library(clusterProfiler))
eKEGG <- enricher(gene = , TERM2GENE = kegg_t2g,
                  universe = ,
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.1, qvalueCutoff = 0.05,
                  minGSSize = 10, maxGSSize = 500)

```

## Differential expression


```r
suppressMessages(library(DESeq2))
ex <- makeExampleDESeqDataSet(m=4)
dds <- DESeq(ex)
#> estimating size factors
#> estimating dispersions
#> gene-wise dispersion estimates
#> mean-dispersion relationship
#> final dispersion estimates
#> fitting model and testing
res <- results(dds, contrast=c("condition","B","A"))
rld <- rlogTransformation(ex, blind=TRUE)
dat <- plotPCA(rld, intgroup=c("condition"),returnData=TRUE)
percentVar <- round(100 * attr(dat,"percentVar"))
suppressMessages(library(ggplot2))
ggplot(dat, aes(PC1, PC2, color=condition, shape=condition)) +
geom_point(size=3) +
xlab(paste0("PC1:",percentVar[1],"% variance")) +
ylab(paste0("PC2:",percentVar[2],"% variance"))
```

![plot of chunk DESeq2](figures/DESeq2-1.png)

```r
ex$condition <- relevel(ex$condition, ref="B")
dds2 <- DESeq(dds)
#> using pre-existing size factors
#> estimating dispersions
#> found already estimated dispersions, replacing these
#> gene-wise dispersion estimates
#> mean-dispersion relationship
#> final dispersion estimates
#> fitting model and testing
res <- results(dds2)
write.csv(as.data.frame(res),file="A_vs_B.csv")
```

## Gene co-expression and network analysis

A simple network is furnished with the `GeneNet` documentation example,


```r
library("GeneNet")
#> Loading required package: corpcor
#> Loading required package: longitudinal
#> Loading required package: fdrtool

# A random network with 20 nodes and 10 percent (=19) edges
true.pcor <- ggm.simulate.pcor(20, 0.1)

# convert to edge list
test.results <- ggm.list.edges(true.pcor)[1:19,]

# Rgraphviz
nlab <- LETTERS[1:20]
gr <- network.make.graph( test.results, nlab)
gr
#> A graphNEL graph with directed edges
#> Number of Nodes = 20 
#> Number of Edges = 38
num.nodes(gr)
#> [1] 20
edge.info(gr)
#> $weight
#>      A~G      A~N      B~G      D~T      D~M      D~G      D~K      E~O 
#>  0.47991 -0.52410  0.41918 -0.17285  0.29350 -0.32985 -0.34516  0.99964 
#>      H~I      H~P      I~J      J~R      K~P      K~R      K~M      L~R 
#>  0.39039  0.72445 -0.58593  0.25108 -0.11564 -0.20274 -0.33219  0.75108 
#>      M~N      M~T      N~S 
#>  0.15757 -0.48960  0.42935 
#> 
#> $dir
#>    A~G    A~N    B~G    D~T    D~M    D~G    D~K    E~O    H~I    H~P    I~J 
#> "none" "none" "none" "none" "none" "none" "none" "none" "none" "none" "none" 
#>    J~R    K~P    K~R    K~M    L~R    M~N    M~T    N~S 
#> "none" "none" "none" "none" "none" "none" "none" "none"
gr2 <- network.make.graph( test.results, nlab, drop.singles=TRUE)
gr2
#> A graphNEL graph with directed edges
#> Number of Nodes = 17 
#> Number of Edges = 38
num.nodes(gr2)
#> [1] 17
edge.info(gr2)
#> $weight
#>      A~G      A~N      B~G      D~T      D~M      D~G      D~K      E~O 
#>  0.47991 -0.52410  0.41918 -0.17285  0.29350 -0.32985 -0.34516  0.99964 
#>      H~I      H~P      I~J      J~R      K~P      K~R      K~M      L~R 
#>  0.39039  0.72445 -0.58593  0.25108 -0.11564 -0.20274 -0.33219  0.75108 
#>      M~N      M~T      N~S 
#>  0.15757 -0.48960  0.42935 
#> 
#> $dir
#>    A~G    A~N    B~G    D~T    D~M    D~G    D~K    E~O    H~I    H~P    I~J 
#> "none" "none" "none" "none" "none" "none" "none" "none" "none" "none" "none" 
#>    J~R    K~P    K~R    K~M    L~R    M~N    M~T    N~S 
#> "none" "none" "none" "none" "none" "none" "none" "none"

# plot network
library("Rgraphviz")
#> Loading required package: graph
#> Loading required package: grid
#> 
#> Attaching package: 'Rgraphviz'
#> The following objects are masked from 'package:IRanges':
#> 
#>     from, to
#> The following objects are masked from 'package:S4Vectors':
#> 
#>     from, to
plot(gr, "fdp")
#> Warning in arrows(tail_from[1], tail_from[2], tail_to[1], tail_to[2], col =
#> edgeColor, : zero-length arrow is of indeterminate angle and so skipped
#> Warning in arrows(head_from[1], head_from[2], head_to[1], head_to[2], col =
#> edgeColor, : zero-length arrow is of indeterminate angle and so skipped
```

![plot of chunk GeneNet](figures/GeneNet-1.png)

```r
plot(gr2, "fdp")
#> Warning in arrows(tail_from[1], tail_from[2], tail_to[1], tail_to[2], col =
#> edgeColor, : zero-length arrow is of indeterminate angle and so skipped

#> Warning in arrows(tail_from[1], tail_from[2], tail_to[1], tail_to[2], col =
#> edgeColor, : zero-length arrow is of indeterminate angle and so skipped
#> Warning in arrows(tail_from[1], tail_from[2], tail_to[1], tail_to[2], col =
#> edgeColor, : zero-length arrow is of indeterminate angle and so skipped

#> Warning in arrows(tail_from[1], tail_from[2], tail_to[1], tail_to[2], col =
#> edgeColor, : zero-length arrow is of indeterminate angle and so skipped
```

![plot of chunk GeneNet](figures/GeneNet-2.png)

and a more involved version

```r
set.seed(123454321)
m <- matrix(runif(2500),50)
r <- cor(m)
g <- as.matrix(r>=0.7)+0
f1 <- Heatmap(r)
f2 <- Heatmap(g)
f <- f1+f2
draw(f) # f2 is somewhat twisted
suppressMessages(library(WGCNA))
pwr <- c(1:10, seq(from = 12, to=30, by=2))
sft <- pickSoftThreshold(dat, powerVector = pwr, verbose = 5)
meg <- moduleEigengenes(t(tpm), color=Colors, softPower=6)
TOM <- TOMsimilarity(adjMatrix)
Tree <- hclust(as.dist(1-TOM), method = "average")
plotDendroAndColors(Tree, colors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
```

## Meta-data

```r
suppressMessages(library(recount3))
hs <- available_projects()
dim(subset(hs,file_source=="gtex"))
annotation_options("human")
blood_rse <- create_rse(subset(hs,project=="BLOOD"))
metadata(blood_rse)
dim(blood_rse)
rowRanges(blood_rse)
colnames(colData(blood_rse))
expand_sra_attributes(blood_rse)
```

---

## A list of Bioconductor/CRAN packages

Package | Description
--------|------------
**Bioconductor** |
AnnotationDbi | AnnotationDb objects and their progeny, methods etc.
org.Hs.eg.db | Conversion of Entrez ID -- gene symbols
EnsDb.Hsapiens.v86 | Exposes an annotation databases generated from Ensembl
ensembldb | Retrieve annotation data from an Ensembl based package
TxDb.Hsapiens.UCSC.hg38.knownGene | Annotation of the human genome
INSPEcT | Quantification of the intronic and exonic gene features and the post-transcriptional regulation analysis
graphite | GRAPH Interaction from pathway topological environment
clusterProfiler | Functional profiles for genes and gene clusters
DESSeq2 | Differential gene expression analysis based on the negative binomial distribution
edgeR | Empirical analysis of digital gene expression
WGCNA | Weighted correlation network analysis
ComplexHeatmap | Make complex heatmaps
recount3 | Interface to uniformly processed RNA-seq data
Pi | Priority index, leveraging genetic evidence to prioritise drug targets at the gene and pathway level
Rgraphiz | Interfaces R with the AT&T graphviz library for plotting R graph objects from the graph package
**CRAN** |
GeneNet | Modeling and Inferring Gene Networks
ggplot2 | Data Visualisations Using the grammar of graphics
pheatmap | results visualisation
plyr | Splitting, applying and combining data
