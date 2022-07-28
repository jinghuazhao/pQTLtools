# CVD1/INF1 protein-HF analysis

Modify `config.yaml` to use `MendelianRandomization` v0.6.0 as with bug-fixing in `workflow/scripts/MR_functions.R`.

This furnishes run on CSD3 without the --use-conda option for `snakemake` as all R packages are indeed available locally.

```bash
module load miniconda3/4.5.1
export mypath=${HOME}/COVID-19/miniconda37
conda activate ${mypath}
snakemake -c all
```

giving `results`/`res_Observational.csv` (observational results) and `res_MR_aggregate.csv` (MR results)

