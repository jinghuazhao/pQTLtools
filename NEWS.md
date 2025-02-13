## pQTLtools 0.4

(2025/2/13)

* Another pass through R CMD check --as-cran
* Revise articles esp. pQTLtools, esse, LocusZoom.js (including lz.html and ACE.html).
* Unify DESCRIPTION/README (finally).
* Add support for MathJax and mermaid.
* Add Caprion coloc scripts.
* Add SummarizedExperiment section to es.Rmd --> esse.Rmd.
* Add csq(), protein_altering_variants (so data/), histo.fyi, aria-label.
* Adopt new hg38 LD blocks for turboman().
* Add spectrum.Rmd.
* Renew call to ldops in novelty_check() in accordance with ieugwasr.
* Replace tests/ with most recent PROC results.

## pQTLtools 0.3

(2024/05/26)

* Activate package logo for the GitHub repository.
* Dedicated .R files for (blocks of) functions.
* Use of scope operator for clarification, e.g., ensembldb::genes().
* Use of OPENGWAS_JWT (.Renviron) from <https://api.opengwas.io/profile/>.
* Suggest IlluminaHumanMethylation450kmanifest, OUTRIDER.
* Fix URL in DESCRIPTION/snakemake.Rmd.

## pQTLtools 0.2

(2024/04/24)

### Milestones

* ***2023.12***. Add LocusZoom.js article.
* ***2023.05***. Test data are fully available with medRxiv post of the SCALLOP paper.
* ***2022.12***. A new package [pQTLdata](https://github.com/jinghuazhao/pQTLdata) is created to hold panel and meta data.
* ***2022.06***. It passes [CRAN](https://cran.r-project.org/) checks with no warning.
* ***2021.02***. A web-driven documentation is now available, <https://jinghuazhao.github.io/pQTLtools/>.

### Package

The information here mirrors the package DESCRIPTION,

* Depends R (>=3.5.0), pQTLdata.
* Import dplyr, gap, ggplot2, Rdpack.
* importFrom utils read.table tail.
* Import from lmm as template, use save(compress='xz').
* Replace ChangeLog with NEWS.md.
* LICENSE.md and README.md.
* Reduce size by `sed -i '/ISSN/d' REFERENCES.bib`.
* Suggest Biobase, BioStrings, GenomeInfoDb, GenomicRanges, IRanges, VariantAnnotation.
* Suggest Roxygen2.
* Suggest biomaRt, bookdown, circlize, cowplot.
* Suggest gap.datasets, gwasvcf, htmlwidgets, httr, ieugwasr.
* Suggest knitr, mclust, meta.
* Suggest openxlsx, plotly, plyr.
* Suggest regione, rgl, rmarkdown, rtracklayer.
* Suggest scatterplot3d, seqminer, stringr.
* GitHub action
* inst/Bioconductor/.
* inst/STRING/change_STRING_colors.py.
* inst/UniProt|PPI/README.md.
* inst/snakemake
* pQTLtools.Rmd, bioconductor.Rmd, es.Rmd, LocusZoom.js.Rmd, snakemake.Rmd and SCALLOP-INF.Rmd articles.
* List publications on pQTLs by Sun et al. (2018) and Suhre et al. (2020).

### Functions

The list is in no particular order,

* peptideMapping(), peptideAssociationPlot().
* make_ExpressionSet(), novelty_check(), qtl_lookup(), turboman(), turboqq().
* Reflow turboman.r/[partial]turboqq.r by formatR::tidy_source().
* pqtlMR(), run_TwoSampleMR().
* run_coloc().
* import_OpenGWAS()
* import_eQTLCatalogue().
* genequries(), regionqueries(), snpqueries().
* uniprot2ids().

## pQTLtools 0.1

* First release
