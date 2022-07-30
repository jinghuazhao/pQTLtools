# CVD1/INF1 protein-HF analysis

Included here is data on the 26 overlapping proteins (no information on IL-4) from the Olink cvd1 & inf1 panels.

Steps to set up the environment are outlined in [notes](notes/README.md), while `MendelianRandomization` v0.6.0 is used together with a bug fix in `workflow/scripts/MR_functions.R`.
The directory `input/` can potentially be built from rules defined in the workflow.


```bash
module load miniconda3/4.5.1
export mypath=${HOME}/COVID-19/miniconda37
conda activate ${mypath}
# A dry run to show details of the workflow.
snakemake -n
# Analysis (no --use-conda option since all R packages are more up-to-date locally)
snakemake -c all
```

which gives `output`/`MR.csv` (MR results) and `Obs.csv` (observational results).
