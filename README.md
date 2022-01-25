# pQTLtools <img src="man/figures/logo.svg" align="right" />

# pQTLtools

## A protein Quantitative Trait Locus toolkit

This seeds collection of data and utilties for (pQTL) analysis. At
this early stage, the repository collects information on a number of
protein panels, linking function for cis/trans classification, 2D
manhattan plots, 3D-plotly plots, forest plots among others availale
from R/gap; query results on genes, regions, and SNPs via
PhenoScanner, adding functionality to check for replication across
platforms and aspects of protein-related analysis such as
pQTL-Mendelian Randomization via TwoSampleMR, linkage through UniProt
IDs to other resources.

## Installation

The latest version of pQTLtools can be installed as usual:

```
install.packages("remotes")
remotes::install_github("jinghuazhao/pQTLtools")
```

Dependencies are detailed in the DECRIPTION file of the package at GitHub.

## February 2021 update 

**A web-driven documentation is now available**

[https://jinghuazhao.github.io/pQTLtools/](https://jinghuazhao.github.io/pQTLtools/)

## A summary of datasets and functions

Objects             |    Description
--------------------|-----------------------------------------
**Datasets**        |    
biomaRt             |    Curated data from biomaRt
caprion             |    Caprion panel
hg19                |    Curated data from Bioconductor
hg19Tables          |    Curated data from UCSC genome browser
inf1                |    Olink/INF panel
Olink_NGS           |    Olink/NGS panels
Olink_qPCR          |    Olink/qPCR panels
SomaLogic160410     |    SomaLogic panel
SomaScanV4.1        |    SomaScan v4.1 panel
st4                 |    ST4 of the INTERVAL SomaLogic paper
st6                 |    ST6 of the INTERVAL SomaLogic paper
st18                |    ST18 of the INTERVAL SomaLogic paper
swath_ms            |    SWATH-MS panel
**eQTL/GWAS**       |
get.prop.below.LLOD  |   Limit of detection analysis
import_eQTLCatalogue |   Import eQTL Catalogue
import_OpenGWAS      |   Import OpenGWAS
make_ExpressionSet   |   A call to ExpressionSet class
run_coloc            |   Colocalisation analysis
**MR analysis**      |
pqtlMR               |   Bidirectional pQTL-MR analysis
run_TwoSampleMR      |   A generic wrapper for TwoSampleMR analysis
**PhenoScanner Utilities** |
genequeries          |   phenoscanner genequeries in batches
regionqueries        |   phenoscanner regionqueries in batches
snpqueries           |   phenoscanner snpqueries in batches
**UniProt API**      |
uniprot2ids          |   UniProt ID to others
