---
title: "Reference data for turboman"
output:
  bookdown::html_document2:
    toc: true
    toc_float: true
    number_sections: true
    self_contained: true
fontsize: 11pt
bibliography: '/rds/project/rds-4o5vpvAowP0/software/R/pQTLtools/REFERENCES.bib'
csl: nature-genetics.csl
pkgdown:
  as_is: true
vignette: >
  %\VignetteIndexEntry{An Overview of pQTLtools}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

## hg19 (build 37) reference data

Details of the reference data are shown as follows,


``` r
options(width=200)
rm(list=ls())
load("~/cambridge-ceu/turboman/turboman_hg19_reference_data.rda")
ls()
```

```
## [1] "ld_block_breaks_pickrell_hg19_eur" "refgene_gene_coordinates_h19"
```

``` r
refgene_gene_coordinates_hg19_eur <- refgene_gene_coordinates_h19
head(ld_block_breaks_pickrell_hg19_eur)
```

```
##   chr   start
## 1   1   10583
## 2   1 1892607
## 3   1 3582736
## 4   1 4380811
## 5   1 5913893
## 6   1 7247335
```

``` r
dim(ld_block_breaks_pickrell_hg19_eur)
```

```
## [1] 1725    2
```

``` r
head(refgene_gene_coordinates_hg19_eur)
```

```
##       chromosome gene_transcription_start gene_transcription_stop  gene_name gene_transcription_midposition
## 1              1                    11873                   14409    DDX11L1                        13141.0
## 53009          1                    17368                   17436  MIR6859-1                        17402.0
## 55877          1                    17368                   17436  MIR6859-3                        17402.0
## 2              1                    14361                   29370     WASH7P                        21865.5
## 45554          1                    30365                   30503 MIR1302-10                        30434.0
## 45166          1                    30365                   30503  MIR1302-9                        30434.0
```

``` r
dim(refgene_gene_coordinates_hg19_eur)
```

```
## [1] 27285     5
```

``` r
save(ld_block_breaks_pickrell_hg19_eur,refgene_gene_coordinates_hg19_eur,
     file="turboman_hg19_reference_data.rda",compress='xz')
```

## liftover

The script is as follows,


``` r
liftover <- function(chr_start_end_snpid)
# chain import currently only handles local, uncompressed file paths
# https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz
# https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz
{
  HPC_WORK <- Sys.getenv("HPC_WORK")
  f <- file.path(HPC_WORK,"bin","hg19ToHg38.over.chain")
  chain <- rtracklayer::import.chain(f)
  names(chain) <- gsub("chr","",names(chain))
  require(GenomicRanges)
  gr <- with(chr_start_end_snpid, GenomicRanges::GRanges(seqnames=chr,IRanges::IRanges(start,end),snpid=snpid))
  seqlevelsStyle(gr) <- "NCBI"
  gr38 <- rtracklayer::liftOver(gr, chain)
  if (all(is.na(gr38$seqnames))) {
    warning("Liftover failed for some SNPs.")
  }
  as.data.frame(gr38)
}

library(dplyr)
library(valr)
ld <- ld_block_breaks_pickrell_hg19_eur %>%
      dplyr::mutate(end=start,snpid=paste(chr,start,end,sep=":")) %>%
      dplyr::select(chr,start,end,snpid)
ld38 <- liftover(ld)
ld_block_breaks_pickrell_hg38_eur <- dplyr::left_join(ld, ld38, by="snpid") %>%
                                     dplyr::transmute(chr,start=start.y) %>%
                                     dplyr::filter(!is.na(start))
head(ld_block_breaks_pickrell_hg38_eur)
```

```
##   chr   start
## 1   1   10583
## 2   1 1961168
## 3   1 3666172
## 4   1 4320751
## 5   1 5853833
## 6   1 7187275
```

``` r
dim(ld_block_breaks_pickrell_hg38_eur)
```

```
## [1] 1723    2
```

## refGene for hg38 (build 38)


``` r
refGene38 <- read.delim("~/tests/turboman/refGene38.tsv") %>%
             setNames(c("chrom","start","end","gene"))
refgene_gene_coordinates_hg38_eur <- valr::bed_merge(dplyr::group_by(refGene38,gene)) %>%
                                     valr::bed_sort() %>%
                                     setNames(c("chromosome", "gene_transcription_start", "gene_transcription_stop", "gene_name")) %>%
                                     dplyr::mutate(chromosome=gsub("chr","",chromosome),
                                                   gene_transcription_midposition=(gene_transcription_start+gene_transcription_stop)/2)
head(refgene_gene_coordinates_hg38_eur)
```

```
## # A tibble: 6 × 5
## # Groups:   gene_name [6]
##   chromosome gene_transcription_start gene_transcription_stop gene_name   gene_transcription_midposition
##   <chr>                         <int>                   <int> <chr>                                <dbl>
## 1 1                             11873                   14409 DDX11L1                             13141 
## 2 1                             14361                   29370 WASH7P                              21866.
## 3 1                             17368                   17436 MIR6859-1                           17402 
## 4 1                             29773                   35418 MIR1302-2HG                         32596.
## 5 1                             30365                   30503 MIR1302-2                           30434 
## 6 1                             34610                   36081 FAM138A                             35346.
```

``` r
dim(refgene_gene_coordinates_hg38_eur)
```

```
## [1] 48463     5
```

## hg38 reference data

LD blocks and refGene are saved into a `.rda` file.

```r
save(ld_block_breaks_pickrell_hg38_eur,refgene_gene_coordinates_hg38_eur,
     file="turboman_hg38_reference_data.rda",compress='xz')
dir(pattern="*rda")
```

## Remarks

The original reference data contains LD blocks as well refGene information and changes have been made so that
1. Only the boundaries of the LD blocks are lifted over from hg19 to hg38.
2. RefGene data would have been updated, so are replaced using external data.
3. For consistency, `refgene_gene_coordinates_h19` is renamed as `refgene_gene_coordinates_hg19_eur`, which would be reflected in the function but a fuzzy prefix `refgene_gene_coordinates_h` would accommodate the original `.rda`.
4. The `elongated` list per gene is merged, such that

chr|start|end|gene
---|-----|---|-----------
1|66533360|66743095|SGIP1
1|66533592|66690375|SGIP1
1|66533592|66743095|SGIP1
1|66534152|66690375|SGIP1
1|66534152|66743095|SGIP1
1|66534152|66729362|SGIP1

becomes `1|66533360|66743095|SGIP1` as composed to `1|66999251|67216822| SGIP1` (hg19).

## External data

This is from <https://raw.githubusercontent.com/jmacdon/LDblocks_GRCh38/master/data/deCODE_EUR_LD_blocks.bed>.


``` r
ld_block_breaks_macdonald_hg38_eur <- read.delim("~/tests/turboman/deCODE_EUR_LD_blocks.bed") %>%
                                      transmute(chr=as.numeric(gsub("chr","",chr)),start)
save(ld_block_breaks_macdonald_hg38_eur,refgene_gene_coordinates_hg38_eur,
     file="turboman_hg38_reference_data.rda",compress='xz')
```
