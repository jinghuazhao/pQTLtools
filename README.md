---
output: github_document
---



<img src="man/figures/logo.svg" align="right" alt="" width="120" />

<!-- badges: start -->
[![pages-build-deployment](https://github.com/jinghuazhao/pQTLtools/actions/workflows/pages/pages-build-deployment/badge.svg)](https://github.com/jinghuazhao/pQTLtools/actions/workflows/pages/pages-build-deployment)
[![R-CMD-check](https://github.com/jinghuazhao/pQTLtools/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jinghuazhao/pQTLtools/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## A protein Quantitative Trait Locus toolkit

It collects data and utilities for pQTL analysis, including 1. Protein GWAS facilities such as Manhattan/QQ/LocusZoom.js plots, novelty/consequence checking; 2. An article linking functions for cis/trans classification, pQTL-gene plot, 2d/3d-plotly plots, forest plots among others available from gap, <https://cran.r-project.org/package=gap>, as well as colocalization, pQTL-Mendelian Randomization via TwoSampleMR. <https://mrcieu.github.io/TwoSampleMR/>; 3. Query on genes, regions, and SNPs via PhenoScanner, <http://www.phenoscanner.medschl.cam.ac.uk/>. 4. Mapping of UniProt IDs to other resources; 5. Showcases of Bioconductor, <https://github.com/bioconductor>, and snakemake, <https://github.com/snakemake>.

## Installation

The latest version of pQTLtools can be installed as usual:

### 1. Install from R

```r
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("jinghuazhao/pQTLtools")
```

### 2. Install from GitHub repository

```bash
git clone https://github.com/jinghuazhao/pQTLtools
R CMD INSTALL pQTLtools
```

Dependencies are detailed in the DECRIPTION file of the package at GitHub.

## A summary of functions

This can be seen from R with

```r
library(help=pQTLtools)
```

or

```r
library(pQTLtools)
?pQTLtools
```
