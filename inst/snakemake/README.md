# CVD1/INF1 protein-HF analysis

File `config.yaml` is modified to use `MendelianRandomization` v0.6.0 and a bug with `workflow/scripts/MR_functions.R` fixed.

The following script furnishes run on CSD3 (without the --use-conda option for `snakemake` as all R packages are indeed available locally).

```bash
module load miniconda3/4.5.1
export mypath=${HOME}/COVID-19/miniconda37
conda activate ${mypath}
snakemake -c all
```

giving `results`/`res_Observational.csv` (observational results) and `res_MR_aggregate.csv` (MR results)

See [work/README.Rmd](work/README.Rmd) for steps to set up the workflow.
