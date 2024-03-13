## pQTLtools 0.2

### Landmarks by date

* ***2023.12***. Add LocusZoom.js vignette
* ***2023.05***. Test data are fully available with medRxiv post of the SCALLOP paper.
* ***2022.12***. A new package [pQTLdata](https://github.com/jinghuazhao/pQTLdata) is created to hold panel and meta data.
* ***2022.06***. It passes [CRAN](https://cran.r-project.org/) checks with no warning.
* ***2021.02***. A web-driven documentation is now available, [https://jinghuazhao.github.io/pQTLtools/](https://jinghuazhao.github.io/pQTLtools/)

### Accumulated changes

* First release.
* Depends R (>=3.5.0), pQTLdata
* Import dplyr, gap, ggplot2, Rdpack
* importFrom utils read.table tail.
* Import from lmm as template, use save(compress='xz').
* Reduce size by `sed -i '/ISSN/d' REFERENCES.bib`
* Suggest htmlwidgets, plotly.
* Suggest bookdown, cowplot, gap.datasets, httr, plyr, rmarkdown, Biobase, stringr.
* Suggest circlize, openxlsx, knitr and add HTML vignette (biomaRt, regioneR).
  with cis/trans-classification/ideogram/mhtplot2d examples.
* Suggest GenomicRanges & IRanges to handle >1MB region in regionqueries with wait=option.
* Replace ChangeLog with NEWS.md and generate vignette/articles with pkgdown.
* Invoke Roxygen2 for documentation.
* Add GitHub action
* Add inst/Bioconductor/.
* Add inst/STRING/change_STRING_colors.py.
* Add inst/UniProt|PPI/README.md.
* Add inst/snakemake
* Add pQTLtools.Rmd, bioconductor.Rmd, es.Rmd, LocusZoom.js.Rmd, snakemake.Rmd and SCALLOP-INF.Rmd articles.
* Add LICENSE.md and README.md.
* Add listed publications on pQTLs by Sun et al. (2018) and Suhre et al. (2020).
* Add run_TwoSampleMR(), make_ExpressionSet(), novelty_check(), qtl_lookup(), turboman, turboqq
* Add pqtlMR() based on TwoSampleMR.
* Add run_coloc() based on coloc.
* Add import_OpenGWAS and therefore suggests gwasvcf, ieugwasr, rtracklayer, VariantAnnotation.
* Add import_eQTLCatalogue from eQTL-Catalogue-resources suggesting seqminer.
* Add genequries, regionqueries, snpqueries.
* Add uniprot2ids() based on UniProt.
