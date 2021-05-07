# pQTLtools 0.1

* ...
* Suggest cowplot, gap.datasets, httr, plyr, rmarkdown, Biobase, stringr
* Suggest circlize, openxlsx, knitr and add HTML vignette (biomaRt, karyoploteR, regioneR)
  with cis/trans-classification/ideogram/mhtplot2d examples
* Suggest GenomicRanges & IRanges to handle >1MB region in regionqueries with wait= option
* Add pQTLtools.Rmd, es.Rmd, gap.Rmd, plogp.Rmd and SCALLOP-INF.Rmd articles
* Add LICENSE.md and README.md
* Add listed publications on pQTLs by Sun et al. (2018) and Suhre et al. (2020)
* Add run_TwoSampleMR(), make_ExpressionSet
* Add pqtlMR() based on TwoSampleMR
* Add run_coloc() based on coloc
* Add import_OpenGWAS and therefore suggests gwasvcf, rtracklayer, VariantAnnotation
* Add import_eQTLCatalogue from eQTL-Catalogue-resources suggesting seqminer
* Add genequries, regionqueries, snpqueries
* Add uniprot2ids() based on UniProt
* Add biomaRt.rda, hg19.rda, hg19Tables.rda, inf1.rda, st4.rda, st6.rda, but drop hgTables.
  Q8NF90 and Q8WWJ7 in inf1.rda were not listed at the UCSC, and replaced with P12034 and P30203 as on UniProt
  Tidy up various options of SomaLogic lookup (panel, box, ST4, ST6)
* Import ggplot2
* Replace ChangeLog with NEWS.md and document with pkgdown.
  Invoke Roxygen2 for documentation
* Import from lmm as template, use save(,compress='xz')
* First release
