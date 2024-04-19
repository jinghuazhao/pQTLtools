## pQTLtools 0.2

### Milestones

* ***2023.12***. Add LocusZoom.js article.
* ***2023.05***. Test data are fully available with medRxiv post of the SCALLOP paper.
* ***2022.12***. A new package [pQTLdata](https://github.com/jinghuazhao/pQTLdata) is created to hold panel and meta data.
* ***2022.06***. It passes [CRAN](https://cran.r-project.org/) checks with no warning.
* ***2021.02***. A web-driven documentation is now available, [https://jinghuazhao.github.io/pQTLtools/](https://jinghuazhao.github.io/pQTLtools/).

### Package

The information here mirrors the package DESCRIPTION,

* Depends R (>=3.5.0), pQTLdata
* Import dplyr, gap, ggplot2, Rdpack
* importFrom utils read.table tail.
* Import from lmm as template, use save(compress='xz').
* Replace ChangeLog with NEWS.md.
* LICENSE.md and README.md.
* Reduce size by `sed -i '/ISSN/d' REFERENCES.bib`
* Suggest htmlwidgets, plotly
* Suggest BioStrings, GenomeInfoDb
* Suggest bookdown, cowplot, gap.datasets, httr, plyr, rmarkdown, Biobase, stringr
* Suggest circlize, openxlsx, knitr and add HTML vignette (biomaRt, regioneR)
* Suggest gwasvcf, ieugwasr, rtracklayer, VariantAnnotation
* Suggest seqminer
* Suggest GenomicRanges, IRanges
* Suggest Roxygen2
* GitHub action
* inst/Bioconductor/.
* inst/STRING/change_STRING_colors.py.
* inst/UniProt|PPI/README.md.
* inst/snakemake
* pQTLtools.Rmd, bioconductor.Rmd, es.Rmd, LocusZoom.js.Rmd, snakemake.Rmd and SCALLOP-INF.Rmd articles.
* List publications on pQTLs by Sun et al. (2018) and Suhre et al. (2020).

### Functions

The list is in no particular order,

* peptideMapping(), peptideAssociationPlot()
* make_ExpressionSet(), novelty_check(), qtl_lookup(), turboman(), turboqq()
* Reflow turboman.r/[partial]turboqq.r by formatR::tidy_source()
* pqtlMR(), run_TwoSampleMR().
* run_coloc().
* import_OpenGWAS()
* import_eQTLCatalogue().
* genequries(), regionqueries(), snpqueries().
* uniprot2ids().

## pQTLtools 0.1

* First release
