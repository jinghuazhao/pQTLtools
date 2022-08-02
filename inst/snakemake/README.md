# Mendelian randomization analysis

[![DOI](https://zenodo.org/badge/429122036.svg)](https://zenodo.org/badge/latestdoi/429122036)
[![Circulation](https://www.ahajournals.org/pb-assets/images/logos/circ-logo-1526571039097.svg)](https://www.ahajournals.org/doi/10.1161/CIRCULATIONAHA.121.056663)

Currently `input/` contains data on the 26 overlapping proteins (no information on IL-4) from the Olink cvd1 & inf1 panels and heart failures.

Steps to set up the environment are outlined in [notes](notes/README.md), while `MendelianRandomization` v0.6.0 is used together with a bug fix in `workflow/r/MR_functions.R`. The workflow has been heavily edited for simplicity and efficiency.


```bash
module load miniconda3/4.5.1
export mypath=${HOME}/COVID-19/miniconda37
source activate ${mypath}
# a dry run (-n).
snakemake --dry-run
# run (-c all on all available cores without --use-conda option as local packages are more up-to-date)
snakemake --cores all
```

The document is knitted with `Rscript -e 'knitr::knit("README.Rmd")'` which also gives `output`/`MR.csv` (MR results) and `Obs.csv` (meta-analysis results based on observational studies).

Some related operations are also ready.

```bash
snakemake --dag | \
dot -Tpdf > dag.pdf
snakemake --rulegraph | \
dot -Tpdf > rulegraph.pdf
snakemake --report report.html
```
