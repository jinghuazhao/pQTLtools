---
title: "Reference data for turboman"
output:
  bookdown::html_document2:
    toc: true
    toc_float: true
    number_sections: true
    self_contained: true
fontsize: 11pt
bibliography: '/rds/project/jmmh2/rds-jmmh2-public_databases/software/R/pQTLtools/REFERENCES.bib'
csl: nature-genetics.csl
pkgdown:
  as_is: true
vignette: >
  %\VignetteIndexEntry{An Overview of pQTLtools}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

## Preparation

The following script has been used to generate the reference data,


```r
options(width=200)
rm(list=ls())
dir()
```

```
## [1] "README.md"                        "README.Rmd"                       "turboman_hg19_reference_data.rda" "urboman_hg38_reference_data.rda"
```

```r
load("~/cambridge-ceu/turboman/turboman_hg19_reference_data.rda")
ls()
```

```
## [1] "ld_block_breaks_pickrell_hg19_eur" "refgene_gene_coordinates_h19"
```

```r
liftover <- function(chr_start_end_snpid)
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
refgene_gene_coordinates_hg19_eur <- refgene_gene_coordinates_h19
save(ld_block_breaks_pickrell_hg19_eur,refgene_gene_coordinates_hg19_eur,file="turboman_hg19_reference_data.rda")
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

```r
ld <- ld_block_breaks_pickrell_hg19_eur %>%
      dplyr::mutate(end=start,snpid=paste(chr,start,end,sep=":")) %>%
      dplyr::select(chr,start,end,snpid)
ld38 <- liftover(ld)
ld_block_breaks_pickrell_hg38_eur <- dplyr::left_join(ld, ld38, by="snpid") %>%
                                     dplyr::transmute(chr,start=start.y) %>%
                                     dplyr::filter(!is.na(start))
refGene38 <- read.delim("~/tests/turboman/refGene38.tsv") %>%
             setNames(c("chrom","start","end","gene"))
refgene_gene_coordinates_hg38_eur <- valr::bed_merge(dplyr::group_by(refGene38,gene)) %>%
                                     valr::bed_sort() %>%
                                     setNames(c("chromosome", "gene_transcription_start", "gene_transcription_stop", "gene_name")) %>%
                                     dplyr::mutate(chromosome=gsub("chr","",chromosome),
                                                   gene_transcription_midposition=(gene_transcription_start+gene_transcription_stop)/2)
save(ld_block_breaks_pickrell_hg38_eur,refgene_gene_coordinates_hg38_eur,file="urboman_hg38_reference_data.rda")
```

Note that the original reference data contains LD blocks as well refGene information and changs have been made so that
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
