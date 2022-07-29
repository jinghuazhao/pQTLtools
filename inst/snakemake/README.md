# CVD1/INF1 protein-HF analysis

Included here is data on a total of 26 in the 27 overlapping proteins from the Olink cvd1 & inf1 panels.

Steps to set up the environment are outlined in [work/README.md](work/README.md). Note `config.yaml` uses `MendelianRandomization` v0.6.0 and a bug in `workflow/scripts/MR_functions.R` is fixed.
The directory `input/` can potentially be built from rules defined in the workflow.


```bash
module load miniconda3/4.5.1
export mypath=${HOME}/COVID-19/miniconda37
conda activate ${mypath}
# A dry run to show details of the workflow.
snakemake -n
# Analysis (no --use-conda option for `snakemake` as all R packages are locally available)
snakemake -c all
```

```
## 
## CommandNotFoundError: Your shell has not been properly configured to use 'conda activate'.
## If your shell is Bash or a Bourne variant, enable conda for the current user with
## 
##     $ echo ". /usr/local/Cluster-Apps/miniconda3/4.5.1/etc/profile.d/conda.sh" >> ~/.bashrc
## 
## or, for all users, enable conda with
## 
##     $ sudo ln -s /usr/local/Cluster-Apps/miniconda3/4.5.1/etc/profile.d/conda.sh /etc/profile.d/conda.sh
## 
## The options above will permanently enable the 'conda' command, but they do NOT
## put conda's base (root) environment on PATH.  To do so, run
## 
##     $ conda activate
## 
## in your terminal, or to put the base environment on PATH permanently, run
## 
##     $ echo "conda activate" >> ~/.bashrc
## 
## Previous to conda 4.4, the recommended way to activate conda was to modify PATH in
## your ~/.bashrc file.  You should manually remove the line that looks like
## 
##     export PATH="/usr/local/Cluster-Apps/miniconda3/4.5.1/bin:$PATH"
## 
## ^^^ The above line should NO LONGER be in your ~/.bashrc file! ^^^
## 
## 
## Building DAG of jobs...
## Nothing to be done (all requested files are present and up to date).
## Building DAG of jobs...
## Nothing to be done (all requested files are present and up to date).
## Complete log: .snakemake/log/2022-07-29T110948.248645.snakemake.log
```

which gives `output`/`Observational.csv` (observational results) and `MR_aggregate.csv` (MR results).
