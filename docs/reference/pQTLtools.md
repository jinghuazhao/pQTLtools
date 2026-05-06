# A reference manual

It collects data and utilities for pQTL analysis, including 1. Protein
GWAS facilities such as Manhattan/QQ/LocusZoom.js plots,
novelty/consequence checking; 2. Articles linking functions for
cis/trans classification, pQTL-gene plot, 2d/3d-plotly plots, forest
plots among others available from 'gap'
(<https://cran.r-project.org/package=gap>) as well as colocalization,
pQTL-Mendelian Randomization via 'TwoSampleMR'
(<https://mrcieu.github.io/TwoSampleMR/>); 3. Query on genes, regions,
and SNPs via 'PhenoScanner' (<https://github.com/phenoscanner/>). 4.
Mapping of UniProt IDs to other resources; 5. Showcases of
'Bioconductor' (<https://github.com/bioconductor>) and 'snakemake'
(<https://github.com/snakemake>).

## Details

Available data and functions are listed in the following table.

|                                    |                                                         |
|------------------------------------|---------------------------------------------------------|
| Objects                            | Description                                             |
| **GWAS**                           |                                                         |
| `turboman`                         | Manhattan plots                                         |
| `turboqq`                          | QQ plots                                                |
|                                    |                                                         |
| **pGWAS**                          |                                                         |
| `csq`                              | Variant consequence                                     |
| `novelty_check`                    | Locus novelty check                                     |
| `peptideAssociationPlot`           | peptide association plot                                |
| `peptideMapping`                   | peptide-to-protein mapping                              |
| `protein_altering_variants`        | Protein Altering Variants (PAVs)                        |
|                                    |                                                         |
| **Expression analysis**            |                                                         |
| `get.prop.below.LLOD`              | Limit of detection analysis                             |
| `import_eQTLCatalogue`             | Import eQTL Catalogue                                   |
| `make_ExpressionSet`               | A call to ExpressionSet class                           |
| `run_coloc`                        | Colocalisation analysis                                 |
|                                    |                                                         |
| **MR analysis**                    |                                                         |
| `import_OpenGWAS`                  | Import OpenGWAS                                         |
| `pqtlMR`                           | Bidirectional pQTL-MR analysis                          |
| `qtl_lookup`                       | QTL lookup                                              |
| `run_TwoSampleMR`                  | A generic wrapper for TwoSampleMR analysis              |
|                                    |                                                         |
| **PhenoScanner Utilities**         |                                                         |
| `genequeries`                      | phenoscanner genequeries in batches                     |
| `regionqueries`                    | phenoscanner regionqueries in batches                   |
| `snpqueries`                       | phenoscanner snpqueries in batches                      |
|                                    |                                                         |
| **UniProt API**                    |                                                         |
| `uniprot2ids`                      | UniProt ID to others                                    |
|                                    |                                                         |
| **Functions in gap**               |                                                         |
| `gap::METAL_forestplot`            | Forest plots from metal analysis                        |
| `gap::ci2ms`                       | Effect size and standard error from confidence interval |
| `gap::cis.vs.trans.classification` | a cis/trans classifier                                  |
| `gap::circos.cis.vs.trans.plot`    | circos plot of cis/trans classification                 |
| `gap::circos.mhtplot`              | circos Manhattan plot with gene annotation              |
| `gap::circos.mhtplot2`             | Another circos Manhattan plot                           |
| `gap::cs`                          | Credible set                                            |
| `gap::get_b_se`                    | Get b and se from AF, n, and z                          |
| `gap::get_pve_se`                  | Get pve and its standard error from n, z                |
| `gap::get_sdy`                     | Get sd(y) from AF, n, b, se                             |
| `gap::mr`                          | Mendelian randomization analysis                        |
| `gap::invnormal`                   | Inverse normal transformation                           |
| `gap::log10p`                      | log10(p) for a standard normal deviate                  |
| `gap::log10pvalue`                 | log10(p) for a P value including its scientific format  |
| `gap::logp`                        | log(p) for a normal deviate                             |
| `gap::mhtplot.trunc`               | Truncated Manhattan plot                                |
| `gap::miamiplot2`                  | Miami plot                                              |
| `gap::mr_forestplot`               | Mendelian Randomization forest plot                     |
| `gap::pvalue`                      | P value for a normal deviate                            |
| `gap::qtlClassifier`               | A QTL cis/trans classifier                              |
| `gap::qtlFinder`                   | Distance-based signal identification                    |
| `gap::qtl2dplot`                   | 2D QTL plot                                             |
| `gap::qtl2dplotly`                 | 2D QTL plotly                                           |
| `gap::qtl3dplotly`                 | 3D QTL plotly                                           |

## Usage

Vignettes on package usage:

-   An Overview of pQTLtools. `vignette("pQTLtools")`.

-   Bioconductor Notes. `vignette("bioconductor")`.

-   LocusZoom.js. `vignette("LocusZoom.js")`.

-   snakemake showcases. `vignette("snakemake")`.

## See also

Useful links:

-   <https://github.com/jinghuazhao/pQTLtools/>

-   <https://jinghuazhao.github.io/pQTLtools/>

-   Report bugs at <https://github.com/jinghuazhao/pQTLtools/issues>

## Author

Jing Hua Zhao in collaboration with other colleagues
