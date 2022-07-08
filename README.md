
<img src="man/figures/logo.svg" align="right" alt="" width="120" />

## A protein Quantitative Trait Locus toolkit

This seeds collection of data and utilties for pQTL analysis. At this
early stage, the repository collects information on a number of protein
panels, linking functions for cis/trans classification, 2D manhattan
plots, 3D-plotly plots, forest plots among others availale from R/gap;
query results on genes, regions, and SNPs via PhenoScanner, adding
functionality to check for replication across platforms and aspects of
protein-related analysis such as pQTL-Mendelian Randomization via
TwoSampleMR, linkage through UniProt IDs to other resources.

Note that some steps have been omitted to avoid uses of external data
and some files in `~/pQTLtools/tests` associated with ongoing projects
have not been made public but it is easily done when ready.

## Installation

The latest version of pQTLtools can be installed as usual:

``` r
install.packages("remotes")
remotes::install_github("jinghuazhao/pQTLtools")
```

Dependencies are detailed in the DECRIPTION file of the package at
GitHub.

## June 2022 update

It passes CRAN checks with no warning.

## February 2021 update

A web-driven documentation is now available.

<https://jinghuazhao.github.io/pQTLtools/>

## A summary of datasets and functions

| Objects                    | Description                                |
|----------------------------|--------------------------------------------|
| **Datasets**               |                                            |
| biomaRt                    | Curated data from biomaRt                  |
| caprion                    | Caprion panel                              |
| hg19                       | Curated data from Bioconductor             |
| hg19Tables                 | Curated data from UCSC genome browser      |
| inf1                       | Olink/INF panel                            |
| Olink\_NGS                 | Olink/NGS panels                           |
| Olink\_qPCR                | Olink/qPCR panels                          |
| SomaLogic160410            | SomaLogic panel                            |
| SomaScanV4.1               | SomaScan v4.1 panel                        |
| st4                        | ST4 of the INTERVAL SomaLogic paper        |
| st6                        | ST6 of the INTERVAL SomaLogic paper        |
| st18                       | ST18 of the INTERVAL SomaLogic paper       |
| swath\_ms                  | SWATH-MS panel                             |
| **eQTL/GWAS**              |                                            |
| get.prop.below.LLOD        | Limit of detection analysis                |
| import\_eQTLCatalogue      | Import eQTL Catalogue                      |
| import\_OpenGWAS           | Import OpenGWAS                            |
| make\_ExpressionSet        | A call to ExpressionSet class              |
| novelty\_check             | Locus novelty check                        |
| run\_coloc                 | Colocalisation analysis                    |
| **MR analysis**            |                                            |
| pqtlMR                     | Bidirectional pQTL-MR analysis             |
| run\_TwoSampleMR           | A generic wrapper for TwoSampleMR analysis |
| **PhenoScanner Utilities** |                                            |
| genequeries                | phenoscanner genequeries in batches        |
| regionqueries              | phenoscanner regionqueries in batches      |
| snpqueries                 | phenoscanner snpqueries in batches         |
| **UniProt API**            |                                            |
| uniprot2ids                | UniProt ID to others                       |

------------------------------------------------------------------------

## Closely related functions in R/gap

| Objects                     | Description                                            |
|-----------------------------|--------------------------------------------------------|
| METAL\_forestplot           | Forest plots from metal analysis                       |
| cis.vs.trans.classification | a cis/trans classifier                                 |
| circos.cis.vs.trans.plot    | circos plot of cis/trans classification                |
| circos.mhtplot              | circos Manhattan plot with gene annotation             |
| circos.mhtplot2             | Another circos Manhattan plot                          |
| cs                          | Credible set                                           |
| get\_b\_se                  | Get b and se from AF, n, and z                         |
| get\_pve\_se                | Get pve and its standard error from n, z               |
| get\_sdy                    | Get sd(y) from AF, n, b, se                            |
| gsmr                        | Mendelian randomization analysis                       |
| invnormal                   | Inverse normal transformation                          |
| log10p                      | log10(p) for a standard normal deviate                 |
| log10pvalue                 | log10(p) for a P value including its scientific format |
| logp                        | log(p) for a normal deviate                            |
| mhtplot.trunc               | Truncated Manhattan plot                               |
| qtlClassifier               | A QTL cis/trans classifier                             |
| qtl2dplot                   | 2D QTL plot                                            |
| qtl2dplotly                 | 2D QTL plotly                                          |
| qtl3dplotly                 | 3D QTL plotly                                          |
