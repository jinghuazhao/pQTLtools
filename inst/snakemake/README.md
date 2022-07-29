# CVD1/INF1 protein-HF analysis

A total of 27 proteins were identified from the Olink/cvd1 panel overlapped with inf1. The following script furnishes run on CSD3 (no --use-conda option for `snakemake` as all R packages are locally available).

## setup


```bash
module load miniconda3/4.5.1
export mypath=${HOME}/COVID-19/miniconda37
conda activate ${mypath}
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
```

## dry-run


```bash
snakemake -n
```

```
## Building DAG of jobs...
## Job stats:
## job             count    min threads    max threads
## ------------  -------  -------------  -------------
## MR_analysis        26              1              1
## aggregate_MR        1              1              1
## all                 1              1              1
## obs_analysis        1              1              1
## total              29              1              1
## 
## 
## [Fri Jul 29 09:40:23 2022]
## rule MR_analysis:
##     input: resources/data_MR.csv, resources/cis200kb_ldrho/CCL20.ld, resources/cis200kb_ldrho/CCL20.snplist
##     output: results/MR/CCL20.csv
##     jobid: 5
##     reason: Missing output files: results/MR/CCL20.csv
##     wildcards: protein=CCL20
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 09:40:23 2022]
## rule MR_analysis:
##     input: resources/data_MR.csv, resources/cis200kb_ldrho/MCP-1.ld, resources/cis200kb_ldrho/MCP-1.snplist
##     output: results/MR/MCP-1.csv
##     jobid: 19
##     reason: Missing output files: results/MR/MCP-1.csv
##     wildcards: protein=MCP-1
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 09:40:23 2022]
## rule MR_analysis:
##     input: resources/data_MR.csv, resources/cis200kb_ldrho/TRAIL.ld, resources/cis200kb_ldrho/TRAIL.snplist
##     output: results/MR/TRAIL.csv
##     jobid: 26
##     reason: Missing output files: results/MR/TRAIL.csv
##     wildcards: protein=TRAIL
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 09:40:23 2022]
## rule MR_analysis:
##     input: resources/data_MR.csv, resources/cis200kb_ldrho/CXCL6.ld, resources/cis200kb_ldrho/CXCL6.snplist
##     output: results/MR/CXCL6.csv
##     jobid: 12
##     reason: Missing output files: results/MR/CXCL6.csv
##     wildcards: protein=CXCL6
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 09:40:23 2022]
## rule MR_analysis:
##     input: resources/data_MR.csv, resources/cis200kb_ldrho/CSF-1.ld, resources/cis200kb_ldrho/CSF-1.snplist
##     output: results/MR/CSF-1.csv
##     jobid: 9
##     reason: Missing output files: results/MR/CSF-1.csv
##     wildcards: protein=CSF-1
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 09:40:23 2022]
## rule MR_analysis:
##     input: resources/data_MR.csv, resources/cis200kb_ldrho/SCF.ld, resources/cis200kb_ldrho/SCF.snplist
##     output: results/MR/SCF.csv
##     jobid: 23
##     reason: Missing output files: results/MR/SCF.csv
##     wildcards: protein=SCF
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 09:40:23 2022]
## rule MR_analysis:
##     input: resources/data_MR.csv, resources/cis200kb_ldrho/IL-18.ld, resources/cis200kb_ldrho/IL-18.snplist
##     output: results/MR/IL-18.csv
##     jobid: 16
##     reason: Missing output files: results/MR/IL-18.csv
##     wildcards: protein=IL-18
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 09:40:23 2022]
## rule MR_analysis:
##     input: resources/data_MR.csv, resources/cis200kb_ldrho/CCL3.ld, resources/cis200kb_ldrho/CCL3.snplist
##     output: results/MR/CCL3.csv
##     jobid: 6
##     reason: Missing output files: results/MR/CCL3.csv
##     wildcards: protein=CCL3
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 09:40:23 2022]
## rule MR_analysis:
##     input: resources/data_MR.csv, resources/cis200kb_ldrho/MMP-1.ld, resources/cis200kb_ldrho/MMP-1.snplist
##     output: results/MR/MMP-1.csv
##     jobid: 20
##     reason: Missing output files: results/MR/MMP-1.csv
##     wildcards: protein=MMP-1
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 09:40:23 2022]
## rule MR_analysis:
##     input: resources/data_MR.csv, resources/cis200kb_ldrho/EN-RAGE.ld, resources/cis200kb_ldrho/EN-RAGE.snplist
##     output: results/MR/EN-RAGE.csv
##     jobid: 13
##     reason: Missing output files: results/MR/EN-RAGE.csv
##     wildcards: protein=EN-RAGE
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 09:40:23 2022]
## rule MR_analysis:
##     input: resources/data_MR.csv, resources/cis200kb_ldrho/TRANCE.ld, resources/cis200kb_ldrho/TRANCE.snplist
##     output: results/MR/TRANCE.csv
##     jobid: 27
##     reason: Missing output files: results/MR/TRANCE.csv
##     wildcards: protein=TRANCE
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 09:40:23 2022]
## rule MR_analysis:
##     input: resources/data_MR.csv, resources/cis200kb_ldrho/CD40.ld, resources/cis200kb_ldrho/CD40.snplist
##     output: results/MR/CD40.csv
##     jobid: 8
##     reason: Missing output files: results/MR/CD40.csv
##     wildcards: protein=CD40
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 09:40:23 2022]
## rule MR_analysis:
##     input: resources/data_MR.csv, resources/cis200kb_ldrho/CX3CL1.ld, resources/cis200kb_ldrho/CX3CL1.snplist
##     output: results/MR/CX3CL1.csv
##     jobid: 10
##     reason: Missing output files: results/MR/CX3CL1.csv
##     wildcards: protein=CX3CL1
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 09:40:23 2022]
## rule MR_analysis:
##     input: resources/data_MR.csv, resources/cis200kb_ldrho/SIRT2.ld, resources/cis200kb_ldrho/SIRT2.snplist
##     output: results/MR/SIRT2.csv
##     jobid: 24
##     reason: Missing output files: results/MR/SIRT2.csv
##     wildcards: protein=SIRT2
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 09:40:23 2022]
## rule MR_analysis:
##     input: resources/data_MR.csv, resources/cis200kb_ldrho/IL-6.ld, resources/cis200kb_ldrho/IL-6.snplist
##     output: results/MR/IL-6.csv
##     jobid: 17
##     reason: Missing output files: results/MR/IL-6.csv
##     wildcards: protein=IL-6
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 09:40:23 2022]
## rule MR_analysis:
##     input: resources/data_MR.csv, resources/cis200kb_ldrho/Beta-NGF.ld, resources/cis200kb_ldrho/Beta-NGF.snplist
##     output: results/MR/Beta-NGF.csv
##     jobid: 3
##     reason: Missing output files: results/MR/Beta-NGF.csv
##     wildcards: protein=Beta-NGF
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 09:40:23 2022]
## rule MR_analysis:
##     input: resources/data_MR.csv, resources/cis200kb_ldrho/CCL4.ld, resources/cis200kb_ldrho/CCL4.snplist
##     output: results/MR/CCL4.csv
##     jobid: 7
##     reason: Missing output files: results/MR/CCL4.csv
##     wildcards: protein=CCL4
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 09:40:23 2022]
## rule MR_analysis:
##     input: resources/data_MR.csv, resources/cis200kb_ldrho/MMP-10.ld, resources/cis200kb_ldrho/MMP-10.snplist
##     output: results/MR/MMP-10.csv
##     jobid: 21
##     reason: Missing output files: results/MR/MMP-10.csv
##     wildcards: protein=MMP-10
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 09:40:23 2022]
## rule obs_analysis:
##     input: resources/data_Observational.csv
##     output: results/Observational.csv
##     jobid: 1
##     reason: Missing output files: results/Observational.csv
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 09:40:23 2022]
## rule MR_analysis:
##     input: resources/data_MR.csv, resources/cis200kb_ldrho/FGF-23.ld, resources/cis200kb_ldrho/FGF-23.snplist
##     output: results/MR/FGF-23.csv
##     jobid: 14
##     reason: Missing output files: results/MR/FGF-23.csv
##     wildcards: protein=FGF-23
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 09:40:23 2022]
## rule MR_analysis:
##     input: resources/data_MR.csv, resources/cis200kb_ldrho/VEGF-A.ld, resources/cis200kb_ldrho/VEGF-A.snplist
##     output: results/MR/VEGF-A.csv
##     jobid: 28
##     reason: Missing output files: results/MR/VEGF-A.csv
##     wildcards: protein=VEGF-A
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 09:40:23 2022]
## rule MR_analysis:
##     input: resources/data_MR.csv, resources/cis200kb_ldrho/CXCL1.ld, resources/cis200kb_ldrho/CXCL1.snplist
##     output: results/MR/CXCL1.csv
##     jobid: 11
##     reason: Missing output files: results/MR/CXCL1.csv
##     wildcards: protein=CXCL1
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 09:40:23 2022]
## rule MR_analysis:
##     input: resources/data_MR.csv, resources/cis200kb_ldrho/TNFSF14.ld, resources/cis200kb_ldrho/TNFSF14.snplist
##     output: results/MR/TNFSF14.csv
##     jobid: 25
##     reason: Missing output files: results/MR/TNFSF14.csv
##     wildcards: protein=TNFSF14
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 09:40:23 2022]
## rule MR_analysis:
##     input: resources/data_MR.csv, resources/cis200kb_ldrho/CASP-8.ld, resources/cis200kb_ldrho/CASP-8.snplist
##     output: results/MR/CASP-8.csv
##     jobid: 4
##     reason: Missing output files: results/MR/CASP-8.csv
##     wildcards: protein=CASP-8
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 09:40:23 2022]
## rule MR_analysis:
##     input: resources/data_MR.csv, resources/cis200kb_ldrho/IL-8.ld, resources/cis200kb_ldrho/IL-8.snplist
##     output: results/MR/IL-8.csv
##     jobid: 18
##     reason: Missing output files: results/MR/IL-8.csv
##     wildcards: protein=IL-8
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 09:40:23 2022]
## rule MR_analysis:
##     input: resources/data_MR.csv, resources/cis200kb_ldrho/OPG.ld, resources/cis200kb_ldrho/OPG.snplist
##     output: results/MR/OPG.csv
##     jobid: 22
##     reason: Missing output files: results/MR/OPG.csv
##     wildcards: protein=OPG
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 09:40:23 2022]
## rule MR_analysis:
##     input: resources/data_MR.csv, resources/cis200kb_ldrho/HGF.ld, resources/cis200kb_ldrho/HGF.snplist
##     output: results/MR/HGF.csv
##     jobid: 15
##     reason: Missing output files: results/MR/HGF.csv
##     wildcards: protein=HGF
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 09:40:23 2022]
## rule aggregate_MR:
##     input: results/MR/Beta-NGF.csv, results/MR/CASP-8.csv, results/MR/CCL20.csv, results/MR/CCL3.csv, results/MR/CCL4.csv, results/MR/CD40.csv, results/MR/CSF-1.csv, results/MR/CX3CL1.csv, results/MR/CXCL1.csv, results/MR/CXCL6.csv, results/MR/EN-RAGE.csv, results/MR/FGF-23.csv, results/MR/HGF.csv, results/MR/IL-18.csv, results/MR/IL-6.csv, results/MR/IL-8.csv, results/MR/MCP-1.csv, results/MR/MMP-1.csv, results/MR/MMP-10.csv, results/MR/OPG.csv, results/MR/SCF.csv, results/MR/SIRT2.csv, results/MR/TNFSF14.csv, results/MR/TRAIL.csv, results/MR/TRANCE.csv, results/MR/VEGF-A.csv
##     output: results/MR_aggregate.csv
##     jobid: 2
##     reason: Missing output files: results/MR_aggregate.csv; Input files updated by another job: results/MR/SIRT2.csv, results/MR/TNFSF14.csv, results/MR/Beta-NGF.csv, results/MR/FGF-23.csv, results/MR/EN-RAGE.csv, results/MR/CCL4.csv, results/MR/CD40.csv, results/MR/IL-6.csv, results/MR/SCF.csv, results/MR/CX3CL1.csv, results/MR/CASP-8.csv, results/MR/MMP-1.csv, results/MR/CCL20.csv, results/MR/TRAIL.csv, results/MR/VEGF-A.csv, results/MR/CCL3.csv, results/MR/IL-8.csv, results/MR/MCP-1.csv, results/MR/MMP-10.csv, results/MR/CXCL1.csv, results/MR/OPG.csv, results/MR/HGF.csv, results/MR/TRANCE.csv, results/MR/IL-18.csv, results/MR/CSF-1.csv, results/MR/CXCL6.csv
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## 
## [Fri Jul 29 09:40:23 2022]
## localrule all:
##     input: results/Observational.csv, results/MR_aggregate.csv
##     jobid: 0
##     reason: Input files updated by another job: results/Observational.csv, results/MR_aggregate.csv
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## Job stats:
## job             count    min threads    max threads
## ------------  -------  -------------  -------------
## MR_analysis        26              1              1
## aggregate_MR        1              1              1
## all                 1              1              1
## obs_analysis        1              1              1
## total              29              1              1
## 
## 
## This was a dry-run (flag -n). The order of jobs does not reflect the order of execution.
```

## run


```bash
snakemake -c all
```

```
## Building DAG of jobs...
## Using shell: /usr/bin/bash
## Provided cores: 32
## Rules claiming more threads will be scaled down.
## Conda environments: ignored
## Job stats:
## job             count    min threads    max threads
## ------------  -------  -------------  -------------
## MR_analysis        26              1              1
## aggregate_MR        1              1              1
## all                 1              1              1
## obs_analysis        1              1              1
## total              29              1              1
## 
## Select jobs to execute...
## 
## [Fri Jul 29 09:40:25 2022]
## rule MR_analysis:
##     input: resources/data_MR.csv, resources/cis200kb_ldrho/CCL20.ld, resources/cis200kb_ldrho/CCL20.snplist
##     output: results/MR/CCL20.csv
##     jobid: 5
##     reason: Missing output files: results/MR/CCL20.csv
##     wildcards: protein=CCL20
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 09:40:25 2022]
## rule MR_analysis:
##     input: resources/data_MR.csv, resources/cis200kb_ldrho/MCP-1.ld, resources/cis200kb_ldrho/MCP-1.snplist
##     output: results/MR/MCP-1.csv
##     jobid: 19
##     reason: Missing output files: results/MR/MCP-1.csv
##     wildcards: protein=MCP-1
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 09:40:25 2022]
## rule MR_analysis:
##     input: resources/data_MR.csv, resources/cis200kb_ldrho/TRAIL.ld, resources/cis200kb_ldrho/TRAIL.snplist
##     output: results/MR/TRAIL.csv
##     jobid: 26
##     reason: Missing output files: results/MR/TRAIL.csv
##     wildcards: protein=TRAIL
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 09:40:25 2022]
## rule MR_analysis:
##     input: resources/data_MR.csv, resources/cis200kb_ldrho/CXCL6.ld, resources/cis200kb_ldrho/CXCL6.snplist
##     output: results/MR/CXCL6.csv
##     jobid: 12
##     reason: Missing output files: results/MR/CXCL6.csv
##     wildcards: protein=CXCL6
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 09:40:25 2022]
## rule MR_analysis:
##     input: resources/data_MR.csv, resources/cis200kb_ldrho/CCL3.ld, resources/cis200kb_ldrho/CCL3.snplist
##     output: results/MR/CCL3.csv
##     jobid: 6
##     reason: Missing output files: results/MR/CCL3.csv
##     wildcards: protein=CCL3
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 09:40:25 2022]
## rule MR_analysis:
##     input: resources/data_MR.csv, resources/cis200kb_ldrho/MMP-1.ld, resources/cis200kb_ldrho/MMP-1.snplist
##     output: results/MR/MMP-1.csv
##     jobid: 20
##     reason: Missing output files: results/MR/MMP-1.csv
##     wildcards: protein=MMP-1
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 09:40:25 2022]
## rule MR_analysis:
##     input: resources/data_MR.csv, resources/cis200kb_ldrho/EN-RAGE.ld, resources/cis200kb_ldrho/EN-RAGE.snplist
##     output: results/MR/EN-RAGE.csv
##     jobid: 13
##     reason: Missing output files: results/MR/EN-RAGE.csv
##     wildcards: protein=EN-RAGE
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 09:40:25 2022]
## rule MR_analysis:
##     input: resources/data_MR.csv, resources/cis200kb_ldrho/TRANCE.ld, resources/cis200kb_ldrho/TRANCE.snplist
##     output: results/MR/TRANCE.csv
##     jobid: 27
##     reason: Missing output files: results/MR/TRANCE.csv
##     wildcards: protein=TRANCE
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 09:40:25 2022]
## rule MR_analysis:
##     input: resources/data_MR.csv, resources/cis200kb_ldrho/CCL4.ld, resources/cis200kb_ldrho/CCL4.snplist
##     output: results/MR/CCL4.csv
##     jobid: 7
##     reason: Missing output files: results/MR/CCL4.csv
##     wildcards: protein=CCL4
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 09:40:25 2022]
## rule MR_analysis:
##     input: resources/data_MR.csv, resources/cis200kb_ldrho/MMP-10.ld, resources/cis200kb_ldrho/MMP-10.snplist
##     output: results/MR/MMP-10.csv
##     jobid: 21
##     reason: Missing output files: results/MR/MMP-10.csv
##     wildcards: protein=MMP-10
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 09:40:25 2022]
## rule obs_analysis:
##     input: resources/data_Observational.csv
##     output: results/Observational.csv
##     jobid: 1
##     reason: Missing output files: results/Observational.csv
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 09:40:25 2022]
## rule MR_analysis:
##     input: resources/data_MR.csv, resources/cis200kb_ldrho/FGF-23.ld, resources/cis200kb_ldrho/FGF-23.snplist
##     output: results/MR/FGF-23.csv
##     jobid: 14
##     reason: Missing output files: results/MR/FGF-23.csv
##     wildcards: protein=FGF-23
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 09:40:25 2022]
## rule MR_analysis:
##     input: resources/data_MR.csv, resources/cis200kb_ldrho/VEGF-A.ld, resources/cis200kb_ldrho/VEGF-A.snplist
##     output: results/MR/VEGF-A.csv
##     jobid: 28
##     reason: Missing output files: results/MR/VEGF-A.csv
##     wildcards: protein=VEGF-A
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 09:40:25 2022]
## rule MR_analysis:
##     input: resources/data_MR.csv, resources/cis200kb_ldrho/OPG.ld, resources/cis200kb_ldrho/OPG.snplist
##     output: results/MR/OPG.csv
##     jobid: 22
##     reason: Missing output files: results/MR/OPG.csv
##     wildcards: protein=OPG
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 09:40:25 2022]
## rule MR_analysis:
##     input: resources/data_MR.csv, resources/cis200kb_ldrho/CD40.ld, resources/cis200kb_ldrho/CD40.snplist
##     output: results/MR/CD40.csv
##     jobid: 8
##     reason: Missing output files: results/MR/CD40.csv
##     wildcards: protein=CD40
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 09:40:25 2022]
## rule MR_analysis:
##     input: resources/data_MR.csv, resources/cis200kb_ldrho/HGF.ld, resources/cis200kb_ldrho/HGF.snplist
##     output: results/MR/HGF.csv
##     jobid: 15
##     reason: Missing output files: results/MR/HGF.csv
##     wildcards: protein=HGF
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 09:40:25 2022]
## rule MR_analysis:
##     input: resources/data_MR.csv, resources/cis200kb_ldrho/CSF-1.ld, resources/cis200kb_ldrho/CSF-1.snplist
##     output: results/MR/CSF-1.csv
##     jobid: 9
##     reason: Missing output files: results/MR/CSF-1.csv
##     wildcards: protein=CSF-1
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 09:40:25 2022]
## rule MR_analysis:
##     input: resources/data_MR.csv, resources/cis200kb_ldrho/SCF.ld, resources/cis200kb_ldrho/SCF.snplist
##     output: results/MR/SCF.csv
##     jobid: 23
##     reason: Missing output files: results/MR/SCF.csv
##     wildcards: protein=SCF
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 09:40:25 2022]
## rule MR_analysis:
##     input: resources/data_MR.csv, resources/cis200kb_ldrho/IL-18.ld, resources/cis200kb_ldrho/IL-18.snplist
##     output: results/MR/IL-18.csv
##     jobid: 16
##     reason: Missing output files: results/MR/IL-18.csv
##     wildcards: protein=IL-18
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 09:40:25 2022]
## rule MR_analysis:
##     input: resources/data_MR.csv, resources/cis200kb_ldrho/CX3CL1.ld, resources/cis200kb_ldrho/CX3CL1.snplist
##     output: results/MR/CX3CL1.csv
##     jobid: 10
##     reason: Missing output files: results/MR/CX3CL1.csv
##     wildcards: protein=CX3CL1
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 09:40:25 2022]
## rule MR_analysis:
##     input: resources/data_MR.csv, resources/cis200kb_ldrho/SIRT2.ld, resources/cis200kb_ldrho/SIRT2.snplist
##     output: results/MR/SIRT2.csv
##     jobid: 24
##     reason: Missing output files: results/MR/SIRT2.csv
##     wildcards: protein=SIRT2
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 09:40:25 2022]
## rule MR_analysis:
##     input: resources/data_MR.csv, resources/cis200kb_ldrho/IL-6.ld, resources/cis200kb_ldrho/IL-6.snplist
##     output: results/MR/IL-6.csv
##     jobid: 17
##     reason: Missing output files: results/MR/IL-6.csv
##     wildcards: protein=IL-6
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 09:40:25 2022]
## rule MR_analysis:
##     input: resources/data_MR.csv, resources/cis200kb_ldrho/Beta-NGF.ld, resources/cis200kb_ldrho/Beta-NGF.snplist
##     output: results/MR/Beta-NGF.csv
##     jobid: 3
##     reason: Missing output files: results/MR/Beta-NGF.csv
##     wildcards: protein=Beta-NGF
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 09:40:25 2022]
## rule MR_analysis:
##     input: resources/data_MR.csv, resources/cis200kb_ldrho/CXCL1.ld, resources/cis200kb_ldrho/CXCL1.snplist
##     output: results/MR/CXCL1.csv
##     jobid: 11
##     reason: Missing output files: results/MR/CXCL1.csv
##     wildcards: protein=CXCL1
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 09:40:25 2022]
## rule MR_analysis:
##     input: resources/data_MR.csv, resources/cis200kb_ldrho/TNFSF14.ld, resources/cis200kb_ldrho/TNFSF14.snplist
##     output: results/MR/TNFSF14.csv
##     jobid: 25
##     reason: Missing output files: results/MR/TNFSF14.csv
##     wildcards: protein=TNFSF14
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 09:40:25 2022]
## rule MR_analysis:
##     input: resources/data_MR.csv, resources/cis200kb_ldrho/CASP-8.ld, resources/cis200kb_ldrho/CASP-8.snplist
##     output: results/MR/CASP-8.csv
##     jobid: 4
##     reason: Missing output files: results/MR/CASP-8.csv
##     wildcards: protein=CASP-8
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 09:40:25 2022]
## rule MR_analysis:
##     input: resources/data_MR.csv, resources/cis200kb_ldrho/IL-8.ld, resources/cis200kb_ldrho/IL-8.snplist
##     output: results/MR/IL-8.csv
##     jobid: 18
##     reason: Missing output files: results/MR/IL-8.csv
##     wildcards: protein=IL-8
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.2 ──
## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.2 ──
## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.2 ──
## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.2 ──
## ✔ ggplot2 3.3.6     ✔ purrr   0.3.4
## ✔ tibble  3.1.8     ✔ dplyr   1.0.9
## ✔ tidyr   1.2.0     ✔ stringr 1.4.0
## ✔ readr   2.1.2     ✔ forcats 0.5.1
## ✔ ggplot2 3.3.6     ✔ purrr   0.3.4
## ✔ tibble  3.1.8     ✔ dplyr   1.0.9
## ✔ tidyr   1.2.0     ✔ stringr 1.4.0
## ✔ readr   2.1.2     ✔ forcats 0.5.1
## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.2 ──
## ✔ ggplot2 3.3.6     ✔ purrr   0.3.4
## ✔ tibble  3.1.8     ✔ dplyr   1.0.9
## ✔ tidyr   1.2.0     ✔ stringr 1.4.0
## ✔ readr   2.1.2     ✔ forcats 0.5.1
## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.2 ──
## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.2 ──
## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.2 ──
## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.2 ──
## ✔ ggplot2 3.3.6     ✔ purrr   0.3.4
## ✔ tibble  3.1.8     ✔ dplyr   1.0.9
## ✔ tidyr   1.2.0     ✔ stringr 1.4.0
## ✔ readr   2.1.2     ✔ forcats 0.5.1
## ✔ ggplot2 3.3.6     ✔ purrr   0.3.4
## ✔ tibble  3.1.8     ✔ dplyr   1.0.9
## ✔ tidyr   1.2.0     ✔ stringr 1.4.0
## ✔ readr   2.1.2     ✔ forcats 0.5.1
## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.2 ──
## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.2 ──
## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.2 ──
## ✔ ggplot2 3.3.6     ✔ purrr   0.3.4
## ✔ tibble  3.1.8     ✔ dplyr   1.0.9
## ✔ tidyr   1.2.0     ✔ stringr 1.4.0
## ✔ readr   2.1.2     ✔ forcats 0.5.1
## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.2 ──
## ✔ ggplot2 3.3.6     ✔ purrr   0.3.4
## ✔ tibble  3.1.8     ✔ dplyr   1.0.9
## ✔ tidyr   1.2.0     ✔ stringr 1.4.0
## ✔ readr   2.1.2     ✔ forcats 0.5.1
## ✔ ggplot2 3.3.6     ✔ purrr   0.3.4
## ✔ tibble  3.1.8     ✔ dplyr   1.0.9
## ✔ tidyr   1.2.0     ✔ stringr 1.4.0
## ✔ readr   2.1.2     ✔ forcats 0.5.1
## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.2 ──
## ✔ ggplot2 3.3.6     ✔ purrr   0.3.4
## ✔ tibble  3.1.8     ✔ dplyr   1.0.9
## ✔ tidyr   1.2.0     ✔ stringr 1.4.0
## ✔ readr   2.1.2     ✔ forcats 0.5.1
## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.2 ──
## ✔ ggplot2 3.3.6     ✔ purrr   0.3.4
## ✔ tibble  3.1.8     ✔ dplyr   1.0.9
## ✔ tidyr   1.2.0     ✔ stringr 1.4.0
## ✔ readr   2.1.2     ✔ forcats 0.5.1
## ✔ ggplot2 3.3.6     ✔ purrr   0.3.4
## ✔ tibble  3.1.8     ✔ dplyr   1.0.9
## ✔ tidyr   1.2.0     ✔ stringr 1.4.0
## ✔ readr   2.1.2     ✔ forcats 0.5.1
## ✔ ggplot2 3.3.6     ✔ purrr   0.3.4
## ✔ tibble  3.1.8     ✔ dplyr   1.0.9
## ✔ tidyr   1.2.0     ✔ stringr 1.4.0
## ✔ readr   2.1.2     ✔ forcats 0.5.1
## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.2 ──
## ✔ ggplot2 3.3.6     ✔ purrr   0.3.4
## ✔ tibble  3.1.8     ✔ dplyr   1.0.9
## ✔ tidyr   1.2.0     ✔ stringr 1.4.0
## ✔ readr   2.1.2     ✔ forcats 0.5.1
## ✔ ggplot2 3.3.6     ✔ purrr   0.3.4
## ✔ tibble  3.1.8     ✔ dplyr   1.0.9
## ✔ tidyr   1.2.0     ✔ stringr 1.4.0
## ✔ readr   2.1.2     ✔ forcats 0.5.1
## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.2 ──
## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.2 ──
## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.2 ──
## ✔ ggplot2 3.3.6     ✔ purrr   0.3.4
## ✔ tibble  3.1.8     ✔ dplyr   1.0.9
## ✔ tidyr   1.2.0     ✔ stringr 1.4.0
## ✔ readr   2.1.2     ✔ forcats 0.5.1
## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.2 ──
## ✔ ggplot2 3.3.6     ✔ purrr   0.3.4
## ✔ tibble  3.1.8     ✔ dplyr   1.0.9
## ✔ tidyr   1.2.0     ✔ stringr 1.4.0
## ✔ readr   2.1.2     ✔ forcats 0.5.1
## ✔ ggplot2 3.3.6     ✔ purrr   0.3.4
## ✔ tibble  3.1.8     ✔ dplyr   1.0.9
## ✔ tidyr   1.2.0     ✔ stringr 1.4.0
## ✔ readr   2.1.2     ✔ forcats 0.5.1
## ✔ ggplot2 3.3.6     ✔ purrr   0.3.4
## ✔ tibble  3.1.8     ✔ dplyr   1.0.9
## ✔ tidyr   1.2.0     ✔ stringr 1.4.0
## ✔ readr   2.1.2     ✔ forcats 0.5.1
## ✔ ggplot2 3.3.6     ✔ purrr   0.3.4
## ✔ tibble  3.1.8     ✔ dplyr   1.0.9
## ✔ tidyr   1.2.0     ✔ stringr 1.4.0
## ✔ readr   2.1.2     ✔ forcats 0.5.1
## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.2 ──
## ✔ ggplot2 3.3.6     ✔ purrr   0.3.4
## ✔ tibble  3.1.8     ✔ dplyr   1.0.9
## ✔ tidyr   1.2.0     ✔ stringr 1.4.0
## ✔ readr   2.1.2     ✔ forcats 0.5.1
## ✔ ggplot2 3.3.6     ✔ purrr   0.3.4
## ✔ tibble  3.1.8     ✔ dplyr   1.0.9
## ✔ tidyr   1.2.0     ✔ stringr 1.4.0
## ✔ readr   2.1.2     ✔ forcats 0.5.1
## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.2 ──
## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.2 ──
## ✔ ggplot2 3.3.6     ✔ purrr   0.3.4
## ✔ tibble  3.1.8     ✔ dplyr   1.0.9
## ✔ tidyr   1.2.0     ✔ stringr 1.4.0
## ✔ readr   2.1.2     ✔ forcats 0.5.1
## ✔ ggplot2 3.3.6     ✔ purrr   0.3.4
## ✔ tibble  3.1.8     ✔ dplyr   1.0.9
## ✔ tidyr   1.2.0     ✔ stringr 1.4.0
## ✔ readr   2.1.2     ✔ forcats 0.5.1
## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.2 ──
## ✔ ggplot2 3.3.6     ✔ purrr   0.3.4
## ✔ tibble  3.1.8     ✔ dplyr   1.0.9
## ✔ tidyr   1.2.0     ✔ stringr 1.4.0
## ✔ readr   2.1.2     ✔ forcats 0.5.1
## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.2 ──
## ✔ ggplot2 3.3.6     ✔ purrr   0.3.4
## ✔ tibble  3.1.8     ✔ dplyr   1.0.9
## ✔ tidyr   1.2.0     ✔ stringr 1.4.0
## ✔ readr   2.1.2     ✔ forcats 0.5.1
## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.2 ──
## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
## ✖ dplyr::filter() masks stats::filter()
## ✖ dplyr::lag()    masks stats::lag()
## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.2 ──
## ✔ ggplot2 3.3.6     ✔ purrr   0.3.4
## ✔ tibble  3.1.8     ✔ dplyr   1.0.9
## ✔ tidyr   1.2.0     ✔ stringr 1.4.0
## ✔ readr   2.1.2     ✔ forcats 0.5.1
## Loading required package: Matrix
## ✔ ggplot2 3.3.6     ✔ purrr   0.3.4
## ✔ tibble  3.1.8     ✔ dplyr   1.0.9
## ✔ tidyr   1.2.0     ✔ stringr 1.4.0
## ✔ readr   2.1.2     ✔ forcats 0.5.1
## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
## ✖ dplyr::between()   masks data.table::between()
## ✖ dplyr::filter()    masks stats::filter()
## ✖ dplyr::first()     masks data.table::first()
## ✖ dplyr::lag()       masks stats::lag()
## ✖ dplyr::last()      masks data.table::last()
## ✖ purrr::transpose() masks data.table::transpose()
## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
## ✖ dplyr::between()   masks data.table::between()
## ✖ dplyr::filter()    masks stats::filter()
## ✖ dplyr::first()     masks data.table::first()
## ✖ dplyr::lag()       masks stats::lag()
## ✖ dplyr::last()      masks data.table::last()
## ✖ purrr::transpose() masks data.table::transpose()
## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
## ✖ dplyr::between()   masks data.table::between()
## ✖ dplyr::filter()    masks stats::filter()
## ✖ dplyr::first()     masks data.table::first()
## ✖ dplyr::lag()       masks stats::lag()
## ✖ dplyr::last()      masks data.table::last()
## ✖ purrr::transpose() masks data.table::transpose()
## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
## ✖ dplyr::between()   masks data.table::between()
## ✖ dplyr::filter()    masks stats::filter()
## ✖ dplyr::first()     masks data.table::first()
## ✖ dplyr::lag()       masks stats::lag()
## ✖ dplyr::last()      masks data.table::last()
## ✖ purrr::transpose() masks data.table::transpose()
## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
## ✖ dplyr::between()   masks data.table::between()
## ✖ dplyr::filter()    masks stats::filter()
## ✖ dplyr::first()     masks data.table::first()
## ✖ dplyr::lag()       masks stats::lag()
## ✖ dplyr::last()      masks data.table::last()
## ✖ purrr::transpose() masks data.table::transpose()
## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
## ✖ dplyr::between()   masks data.table::between()
## ✖ dplyr::filter()    masks stats::filter()
## ✖ dplyr::first()     masks data.table::first()
## ✖ dplyr::lag()       masks stats::lag()
## ✖ dplyr::last()      masks data.table::last()
## ✖ purrr::transpose() masks data.table::transpose()
## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
## ✖ dplyr::between()   masks data.table::between()
## ✖ dplyr::filter()    masks stats::filter()
## ✖ dplyr::first()     masks data.table::first()
## ✖ dplyr::lag()       masks stats::lag()
## ✖ dplyr::last()      masks data.table::last()
## ✖ purrr::transpose() masks data.table::transpose()
## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
## ✖ dplyr::between()   masks data.table::between()
## ✖ dplyr::filter()    masks stats::filter()
## ✖ dplyr::first()     masks data.table::first()
## ✖ dplyr::lag()       masks stats::lag()
## ✖ dplyr::last()      masks data.table::last()
## ✖ purrr::transpose() masks data.table::transpose()
## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
## ✖ dplyr::between()   masks data.table::between()
## ✖ dplyr::filter()    masks stats::filter()
## ✖ dplyr::first()     masks data.table::first()
## ✖ dplyr::lag()       masks stats::lag()
## ✖ dplyr::last()      masks data.table::last()
## ✖ purrr::transpose() masks data.table::transpose()
## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
## ✖ dplyr::between()   masks data.table::between()
## ✖ dplyr::filter()    masks stats::filter()
## ✖ dplyr::first()     masks data.table::first()
## ✖ dplyr::lag()       masks stats::lag()
## ✖ dplyr::last()      masks data.table::last()
## ✖ purrr::transpose() masks data.table::transpose()
## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
## ✖ dplyr::between()   masks data.table::between()
## ✖ dplyr::filter()    masks stats::filter()
## ✖ dplyr::first()     masks data.table::first()
## ✖ dplyr::lag()       masks stats::lag()
## ✖ dplyr::last()      masks data.table::last()
## ✖ purrr::transpose() masks data.table::transpose()
## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
## ✖ dplyr::between()   masks data.table::between()
## ✖ dplyr::filter()    masks stats::filter()
## ✖ dplyr::first()     masks data.table::first()
## ✖ dplyr::lag()       masks stats::lag()
## ✖ dplyr::last()      masks data.table::last()
## ✖ purrr::transpose() masks data.table::transpose()
## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
## ✖ dplyr::between()   masks data.table::between()
## ✖ dplyr::filter()    masks stats::filter()
## ✖ dplyr::first()     masks data.table::first()
## ✖ dplyr::lag()       masks stats::lag()
## ✖ dplyr::last()      masks data.table::last()
## ✖ purrr::transpose() masks data.table::transpose()
## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
## ✖ dplyr::between()   masks data.table::between()
## ✖ dplyr::filter()    masks stats::filter()
## ✖ dplyr::first()     masks data.table::first()
## ✖ dplyr::lag()       masks stats::lag()
## ✖ dplyr::last()      masks data.table::last()
## ✖ purrr::transpose() masks data.table::transpose()
## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
## ✖ dplyr::between()   masks data.table::between()
## ✖ dplyr::filter()    masks stats::filter()
## ✖ dplyr::first()     masks data.table::first()
## ✖ dplyr::lag()       masks stats::lag()
## ✖ dplyr::last()      masks data.table::last()
## ✖ purrr::transpose() masks data.table::transpose()
## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
## ✖ dplyr::between()   masks data.table::between()
## ✖ dplyr::filter()    masks stats::filter()
## ✖ dplyr::first()     masks data.table::first()
## ✖ dplyr::lag()       masks stats::lag()
## ✖ dplyr::last()      masks data.table::last()
## ✖ purrr::transpose() masks data.table::transpose()
## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
## ✖ dplyr::between()   masks data.table::between()
## ✖ dplyr::filter()    masks stats::filter()
## ✖ dplyr::first()     masks data.table::first()
## ✖ dplyr::lag()       masks stats::lag()
## ✖ dplyr::last()      masks data.table::last()
## ✖ purrr::transpose() masks data.table::transpose()
## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
## ✖ dplyr::between()   masks data.table::between()
## ✖ dplyr::filter()    masks stats::filter()
## ✖ dplyr::first()     masks data.table::first()
## ✖ dplyr::lag()       masks stats::lag()
## ✖ dplyr::last()      masks data.table::last()
## ✖ purrr::transpose() masks data.table::transpose()
## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
## ✖ dplyr::between()   masks data.table::between()
## ✖ dplyr::filter()    masks stats::filter()
## ✖ dplyr::first()     masks data.table::first()
## ✖ dplyr::lag()       masks stats::lag()
## ✖ dplyr::last()      masks data.table::last()
## ✖ purrr::transpose() masks data.table::transpose()
## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
## ✖ dplyr::between()   masks data.table::between()
## ✖ dplyr::filter()    masks stats::filter()
## ✖ dplyr::first()     masks data.table::first()
## ✖ dplyr::lag()       masks stats::lag()
## ✖ dplyr::last()      masks data.table::last()
## ✖ purrr::transpose() masks data.table::transpose()
## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
## ✖ dplyr::between()   masks data.table::between()
## ✖ dplyr::filter()    masks stats::filter()
## ✖ dplyr::first()     masks data.table::first()
## ✖ dplyr::lag()       masks stats::lag()
## ✖ dplyr::last()      masks data.table::last()
## ✖ purrr::transpose() masks data.table::transpose()
## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
## ✖ dplyr::between()   masks data.table::between()
## ✖ dplyr::filter()    masks stats::filter()
## ✖ dplyr::first()     masks data.table::first()
## ✖ dplyr::lag()       masks stats::lag()
## ✖ dplyr::last()      masks data.table::last()
## ✖ purrr::transpose() masks data.table::transpose()
## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
## ✖ dplyr::between()   masks data.table::between()
## ✖ dplyr::filter()    masks stats::filter()
## ✖ dplyr::first()     masks data.table::first()
## ✖ dplyr::lag()       masks stats::lag()
## ✖ dplyr::last()      masks data.table::last()
## ✖ purrr::transpose() masks data.table::transpose()
## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
## ✖ dplyr::between()   masks data.table::between()
## ✖ dplyr::filter()    masks stats::filter()
## ✖ dplyr::first()     masks data.table::first()
## ✖ dplyr::lag()       masks stats::lag()
## ✖ dplyr::last()      masks data.table::last()
## ✖ purrr::transpose() masks data.table::transpose()
## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
## ✖ dplyr::between()   masks data.table::between()
## ✖ dplyr::filter()    masks stats::filter()
## ✖ dplyr::first()     masks data.table::first()
## ✖ dplyr::lag()       masks stats::lag()
## ✖ dplyr::last()      masks data.table::last()
## ✖ purrr::transpose() masks data.table::transpose()
## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
## ✖ dplyr::between()   masks data.table::between()
## ✖ dplyr::filter()    masks stats::filter()
## ✖ dplyr::first()     masks data.table::first()
## ✖ dplyr::lag()       masks stats::lag()
## ✖ dplyr::last()      masks data.table::last()
## ✖ purrr::transpose() masks data.table::transpose()
## 
## Attaching package: ‘Matrix’
## 
## The following objects are masked from ‘package:tidyr’:
## 
##     expand, pack, unpack
## 
## Loading required package: metadat
## 
## Loading the 'metafor' package (version 3.4-0). For an
## introduction to the package please type: help(metafor)
## 
## Rows: 128 Columns: 5
## ── Column specification ────────────────────────────────────────────────────────
## Delimiter: ","
## chr (2): Study, Protein
## dbl (3): beta, se, P
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
## Warning messages:
## 1: Problem while computing `meta = map(data, run_meta)`.
## ℹ Studies with NAs omitted from model fitting.
## ℹ The warning occurred in group 1: Protein = "Beta-NGF". 
## 2: Problem while computing `meta = map(data, run_meta)`.
## ℹ Studies with NAs omitted from model fitting.
## ℹ The warning occurred in group 17: Protein = "IL-4". 
## 3: Problem while computing `meta = map(data, run_meta)`.
## ℹ Studies with NAs omitted from model fitting.
## ℹ The warning occurred in group 27: Protein = "SIRT2". 
## 4: `cols` is now required when using unnest().
## Please use `cols = c(meta)` 
## [Fri Jul 29 09:40:30 2022]
## Finished job 1.
## 1 of 29 steps (3%) done
## Method requires data on >2 variants.Method requires data on >2 variants.Method requires data on >2 variants.Method requires data on >2 variants.Method requires data on >2 variants.Method requires data on >2 variants.Method requires data on >2 variants.Method requires data on >2 variants.Method requires data on >2 variants.Method requires data on >2 variants.Method requires data on >2 variants.Method requires data on >2 variants.Method requires data on >2 variants.Method requires data on >2 variants.Method requires data on >2 variants.Method requires data on >2 variants.Method requires data on >2 variants.Method requires data on >2 variants.Method requires data on >2 variants.Method requires data on >2 variants.Method requires data on >2 variants.Method requires data on >2 variants.Method requires data on >2 variants.Method requires data on >2 variants.Method requires data on >2 variants.Method requires data on >2 variants.Method requires data on >2 variants.[Fri Jul 29 09:40:31 2022]
## Finished job 23.
## 2 of 29 steps (7%) done
## Method requires data on >2 variants.Method requires data on >2 variants.Method requires data on >2 variants.Method requires data on >2 variants.Method requires data on >2 variants.Method requires data on >2 variants.Method requires data on >2 variants.[Fri Jul 29 09:40:31 2022]
## Finished job 4.
## 3 of 29 steps (10%) done
## [Fri Jul 29 09:40:31 2022]
## Finished job 18.
## 4 of 29 steps (14%) done
## [Fri Jul 29 09:40:31 2022]
## Finished job 16.
## 5 of 29 steps (17%) done
## [Fri Jul 29 09:40:31 2022]
## Finished job 12.
## 6 of 29 steps (21%) done
## [Fri Jul 29 09:40:31 2022]
## Finished job 7.
## 7 of 29 steps (24%) done
## [Fri Jul 29 09:40:31 2022]
## Finished job 11.
## 8 of 29 steps (28%) done
## Method requires data on >2 variants.Warning messages:
## 1: Problem while computing `MR = map(data, run_MR_all, ldrho)`.
## ℹ NaNs produced
## ℹ The warning occurred in group 23: r2_thresh = 0.4, P_thresh = 0.01, Protein =
##   "TRANCE". 
## 2: Problem while computing `MR = map(data, run_MR_all, ldrho)`.
## ℹ NaNs produced
## ℹ The warning occurred in group 23: r2_thresh = 0.4, P_thresh = 0.01, Protein =
##   "TRANCE". 
## [Fri Jul 29 09:40:32 2022]
## Finished job 17.
## 9 of 29 steps (31%) done
## Warning messages:
## 1: Problem while computing `MR = map(data, run_MR_all, ldrho)`.
## ℹ NaNs produced
## ℹ The warning occurred in group 30: r2_thresh = 0.6, P_thresh = 1, Protein =
##   "CCL3". 
## 2: Problem while computing `MR = map(data, run_MR_all, ldrho)`.
## ℹ NaNs produced
## ℹ The warning occurred in group 30: r2_thresh = 0.6, P_thresh = 1, Protein =
##   "CCL3". 
## Warning messages:
## 1: Problem while computing `MR = map(data, run_MR_all, ldrho)`.
## ℹ NaNs produced
## ℹ The warning occurred in group 30: r2_thresh = 0.6, P_thresh = 1, Protein =
##   "EN-RAGE". 
## 2: Problem while computing `MR = map(data, run_MR_all, ldrho)`.
## ℹ NaNs produced
## ℹ The warning occurred in group 30: r2_thresh = 0.6, P_thresh = 1, Protein =
##   "EN-RAGE". 
## 3: Problem while computing `MR = map(data, run_MR_all, ldrho)`.
## ℹ NaNs produced
## ℹ The warning occurred in group 30: r2_thresh = 0.6, P_thresh = 1, Protein =
##   "EN-RAGE". 
## 4: Problem while computing `MR = map(data, run_MR_all, ldrho)`.
## ℹ NaNs produced
## ℹ The warning occurred in group 30: r2_thresh = 0.6, P_thresh = 1, Protein =
##   "EN-RAGE". 
## 5: Problem while computing `MR = map(data, run_MR_all, ldrho)`.
## ℹ NaNs produced
## ℹ The warning occurred in group 30: r2_thresh = 0.6, P_thresh = 1, Protein =
##   "EN-RAGE". 
## [Fri Jul 29 09:40:32 2022]
## Finished job 27.
## 10 of 29 steps (34%) done
## [Fri Jul 29 09:40:32 2022]
## Finished job 13.
## 11 of 29 steps (38%) done
## [Fri Jul 29 09:40:32 2022]
## Finished job 6.
## 12 of 29 steps (41%) done
## Warning messages:
## 1: Problem while computing `MR = map(data, run_MR_all, ldrho)`.
## ℹ NaNs produced
## ℹ The warning occurred in group 30: r2_thresh = 0.6, P_thresh = 1, Protein =
##   "OPG". 
## 2: Problem while computing `MR = map(data, run_MR_all, ldrho)`.
## ℹ NaNs produced
## ℹ The warning occurred in group 30: r2_thresh = 0.6, P_thresh = 1, Protein =
##   "OPG". 
## 3: Problem while computing `MR = map(data, run_MR_all, ldrho)`.
## ℹ NaNs produced
## ℹ The warning occurred in group 30: r2_thresh = 0.6, P_thresh = 1, Protein =
##   "OPG". 
## 4: Problem while computing `MR = map(data, run_MR_all, ldrho)`.
## ℹ NaNs produced
## ℹ The warning occurred in group 30: r2_thresh = 0.6, P_thresh = 1, Protein =
##   "OPG". 
## Warning messages:
## 1: Problem while computing `MR = map(data, run_MR_all, ldrho)`.
## ℹ NaNs produced
## ℹ The warning occurred in group 30: r2_thresh = 0.6, P_thresh = 1, Protein =
##   "MCP-1". 
## 2: Problem while computing `MR = map(data, run_MR_all, ldrho)`.
## ℹ NaNs produced
## ℹ The warning occurred in group 30: r2_thresh = 0.6, P_thresh = 1, Protein =
##   "MCP-1". 
## 3: Problem while computing `MR = map(data, run_MR_all, ldrho)`.
## ℹ NaNs produced
## ℹ The warning occurred in group 30: r2_thresh = 0.6, P_thresh = 1, Protein =
##   "MCP-1". 
## Warning message:
## Problem while computing `MR = map(data, run_MR_all, ldrho)`.
## ℹ NaNs produced
## ℹ The warning occurred in group 30: r2_thresh = 0.6, P_thresh = 1, Protein =
##   "HGF". 
## [Fri Jul 29 09:40:32 2022]
## Finished job 22.
## 13 of 29 steps (45%) done
## [Fri Jul 29 09:40:32 2022]
## Finished job 19.
## 14 of 29 steps (48%) done
## [Fri Jul 29 09:40:32 2022]
## Finished job 15.
## 15 of 29 steps (52%) done
## Warning messages:
## 1: Problem while computing `MR = map(data, run_MR_all, ldrho)`.
## ℹ NaNs produced
## ℹ The warning occurred in group 25: r2_thresh = 0.6, P_thresh = 1, Protein =
##   "SIRT2". 
## 2: Problem while computing `MR = map(data, run_MR_all, ldrho)`.
## ℹ NaNs produced
## ℹ The warning occurred in group 25: r2_thresh = 0.6, P_thresh = 1, Protein =
##   "SIRT2". 
## 3: Problem while computing `MR = map(data, run_MR_all, ldrho)`.
## ℹ NaNs produced
## ℹ The warning occurred in group 25: r2_thresh = 0.6, P_thresh = 1, Protein =
##   "SIRT2". 
## 4: Problem while computing `MR = map(data, run_MR_all, ldrho)`.
## ℹ NaNs produced
## ℹ The warning occurred in group 25: r2_thresh = 0.6, P_thresh = 1, Protein =
##   "SIRT2". 
## [Fri Jul 29 09:40:32 2022]
## Finished job 24.
## 16 of 29 steps (55%) done
## Warning message:
## Problem while computing `MR = map(data, run_MR_all, ldrho)`.
## ℹ NaNs produced
## ℹ The warning occurred in group 25: r2_thresh = 0.6, P_thresh = 5e-08, Protein
##   = "CX3CL1". 
## [Fri Jul 29 09:40:32 2022]
## Finished job 10.
## 17 of 29 steps (59%) done
## [Fri Jul 29 09:40:33 2022]
## Finished job 9.
## 18 of 29 steps (62%) done
## [Fri Jul 29 09:40:33 2022]
## Finished job 8.
## 19 of 29 steps (66%) done
## [Fri Jul 29 09:40:33 2022]
## Finished job 28.
## 20 of 29 steps (69%) done
## [Fri Jul 29 09:40:33 2022]
## Finished job 3.
## 21 of 29 steps (72%) done
## Warning messages:
## 1: Problem while computing `MR = map(data, run_MR_all, ldrho)`.
## ℹ NaNs produced
## ℹ The warning occurred in group 30: r2_thresh = 0.6, P_thresh = 1, Protein =
##   "FGF-23". 
## 2: Problem while computing `MR = map(data, run_MR_all, ldrho)`.
## ℹ NaNs produced
## ℹ The warning occurred in group 30: r2_thresh = 0.6, P_thresh = 1, Protein =
##   "FGF-23". 
## [Fri Jul 29 09:40:33 2022]
## Finished job 14.
## 22 of 29 steps (76%) done
## Warning messages:
## 1: Problem while computing `MR = map(data, run_MR_all, ldrho)`.
## ℹ NaNs produced
## ℹ The warning occurred in group 24: r2_thresh = 0.4, P_thresh = 1, Protein =
##   "CCL20". 
## 2: Problem while computing `MR = map(data, run_MR_all, ldrho)`.
## ℹ NaNs produced
## ℹ The warning occurred in group 30: r2_thresh = 0.6, P_thresh = 1, Protein =
##   "CCL20". 
## [Fri Jul 29 09:40:34 2022]
## Finished job 5.
## 23 of 29 steps (79%) done
## [Fri Jul 29 09:40:34 2022]
## Finished job 20.
## 24 of 29 steps (83%) done
## Warning messages:
## 1: Problem while computing `MR = map(data, run_MR_all, ldrho)`.
## ℹ NaNs produced
## ℹ The warning occurred in group 30: r2_thresh = 0.6, P_thresh = 1, Protein =
##   "TRAIL". 
## 2: Problem while computing `MR = map(data, run_MR_all, ldrho)`.
## ℹ NaNs produced
## ℹ The warning occurred in group 30: r2_thresh = 0.6, P_thresh = 1, Protein =
##   "TRAIL". 
## 3: Problem while computing `MR = map(data, run_MR_all, ldrho)`.
## ℹ NaNs produced
## ℹ The warning occurred in group 30: r2_thresh = 0.6, P_thresh = 1, Protein =
##   "TRAIL". 
## [Fri Jul 29 09:40:35 2022]
## Finished job 26.
## 25 of 29 steps (86%) done
## Warning message:
## Problem while computing `MR = map(data, run_MR_all, ldrho)`.
## ℹ NaNs produced
## ℹ The warning occurred in group 24: r2_thresh = 0.4, P_thresh = 1, Protein =
##   "MMP-10". 
## [Fri Jul 29 09:40:35 2022]
## Finished job 21.
## 26 of 29 steps (90%) done
## [Fri Jul 29 09:40:36 2022]
## Finished job 25.
## 27 of 29 steps (93%) done
## Select jobs to execute...
## 
## [Fri Jul 29 09:40:36 2022]
## rule aggregate_MR:
##     input: results/MR/Beta-NGF.csv, results/MR/CASP-8.csv, results/MR/CCL20.csv, results/MR/CCL3.csv, results/MR/CCL4.csv, results/MR/CD40.csv, results/MR/CSF-1.csv, results/MR/CX3CL1.csv, results/MR/CXCL1.csv, results/MR/CXCL6.csv, results/MR/EN-RAGE.csv, results/MR/FGF-23.csv, results/MR/HGF.csv, results/MR/IL-18.csv, results/MR/IL-6.csv, results/MR/IL-8.csv, results/MR/MCP-1.csv, results/MR/MMP-1.csv, results/MR/MMP-10.csv, results/MR/OPG.csv, results/MR/SCF.csv, results/MR/SIRT2.csv, results/MR/TNFSF14.csv, results/MR/TRAIL.csv, results/MR/TRANCE.csv, results/MR/VEGF-A.csv
##     output: results/MR_aggregate.csv
##     jobid: 2
##     reason: Missing output files: results/MR_aggregate.csv; Input files updated by another job: results/MR/VEGF-A.csv, results/MR/SIRT2.csv, results/MR/SCF.csv, results/MR/CD40.csv, results/MR/CCL3.csv, results/MR/CXCL1.csv, results/MR/FGF-23.csv, results/MR/CXCL6.csv, results/MR/MMP-10.csv, results/MR/CCL4.csv, results/MR/IL-8.csv, results/MR/CSF-1.csv, results/MR/MMP-1.csv, results/MR/IL-18.csv, results/MR/HGF.csv, results/MR/IL-6.csv, results/MR/CASP-8.csv, results/MR/Beta-NGF.csv, results/MR/CCL20.csv, results/MR/CX3CL1.csv, results/MR/MCP-1.csv, results/MR/TRAIL.csv, results/MR/OPG.csv, results/MR/TNFSF14.csv, results/MR/TRANCE.csv, results/MR/EN-RAGE.csv
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 09:40:36 2022]
## Finished job 2.
## 28 of 29 steps (97%) done
## Select jobs to execute...
## 
## [Fri Jul 29 09:40:36 2022]
## localrule all:
##     input: results/Observational.csv, results/MR_aggregate.csv
##     jobid: 0
##     reason: Input files updated by another job: results/MR_aggregate.csv, results/Observational.csv
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 09:40:36 2022]
## Finished job 0.
## 29 of 29 steps (100%) done
## Complete log: .snakemake/log/2022-07-29T094024.525951.snakemake.log
```

which gives `results`/`Observational.csv` (observational results) and `MR_aggregate.csv` (MR results).

## Additional information

See [work/README.Rmd](work/README.Rmd) for steps to set up the environment, and `config.yaml` is modified to use `MendelianRandomization` v0.6.0 with a bug in `workflow/scripts/MR_functions.R` being fixed.
