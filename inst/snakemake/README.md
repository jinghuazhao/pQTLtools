# CVD1/INF1 protein-HF analysis

Included here is data on the 26 overlapping proteins (no information on IL-4) from the Olink cvd1 & inf1 panels.

Steps to set up the environment are outlined in [notes/README.md](notes/README.md), while `MendelianRandomization` v0.6.0 is used together with a bug fix in `workflow/scripts/MR_functions.R`.
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
## [Fri Jul 29 13:04:22 2022]
## rule MR_analysis:
##     input: input/MR.csv, input/ld/CCL20.ld, input/ld/CCL20.snplist
##     output: output/MR/CCL20.csv
##     jobid: 5
##     reason: Missing output files: output/MR/CCL20.csv
##     wildcards: protein=CCL20
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 13:04:22 2022]
## rule MR_analysis:
##     input: input/MR.csv, input/ld/MCP-1.ld, input/ld/MCP-1.snplist
##     output: output/MR/MCP-1.csv
##     jobid: 19
##     reason: Missing output files: output/MR/MCP-1.csv
##     wildcards: protein=MCP-1
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 13:04:22 2022]
## rule MR_analysis:
##     input: input/MR.csv, input/ld/TRAIL.ld, input/ld/TRAIL.snplist
##     output: output/MR/TRAIL.csv
##     jobid: 26
##     reason: Missing output files: output/MR/TRAIL.csv
##     wildcards: protein=TRAIL
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 13:04:22 2022]
## rule MR_analysis:
##     input: input/MR.csv, input/ld/CXCL6.ld, input/ld/CXCL6.snplist
##     output: output/MR/CXCL6.csv
##     jobid: 12
##     reason: Missing output files: output/MR/CXCL6.csv
##     wildcards: protein=CXCL6
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 13:04:22 2022]
## rule MR_analysis:
##     input: input/MR.csv, input/ld/CSF-1.ld, input/ld/CSF-1.snplist
##     output: output/MR/CSF-1.csv
##     jobid: 9
##     reason: Missing output files: output/MR/CSF-1.csv
##     wildcards: protein=CSF-1
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 13:04:22 2022]
## rule MR_analysis:
##     input: input/MR.csv, input/ld/SCF.ld, input/ld/SCF.snplist
##     output: output/MR/SCF.csv
##     jobid: 23
##     reason: Missing output files: output/MR/SCF.csv
##     wildcards: protein=SCF
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 13:04:22 2022]
## rule MR_analysis:
##     input: input/MR.csv, input/ld/IL-18.ld, input/ld/IL-18.snplist
##     output: output/MR/IL-18.csv
##     jobid: 16
##     reason: Missing output files: output/MR/IL-18.csv
##     wildcards: protein=IL-18
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 13:04:22 2022]
## rule MR_analysis:
##     input: input/MR.csv, input/ld/CCL3.ld, input/ld/CCL3.snplist
##     output: output/MR/CCL3.csv
##     jobid: 6
##     reason: Missing output files: output/MR/CCL3.csv
##     wildcards: protein=CCL3
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 13:04:22 2022]
## rule MR_analysis:
##     input: input/MR.csv, input/ld/MMP-1.ld, input/ld/MMP-1.snplist
##     output: output/MR/MMP-1.csv
##     jobid: 20
##     reason: Missing output files: output/MR/MMP-1.csv
##     wildcards: protein=MMP-1
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 13:04:22 2022]
## rule MR_analysis:
##     input: input/MR.csv, input/ld/TRANCE.ld, input/ld/TRANCE.snplist
##     output: output/MR/TRANCE.csv
##     jobid: 27
##     reason: Missing output files: output/MR/TRANCE.csv
##     wildcards: protein=TRANCE
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 13:04:22 2022]
## rule MR_analysis:
##     input: input/MR.csv, input/ld/EN-RAGE.ld, input/ld/EN-RAGE.snplist
##     output: output/MR/EN-RAGE.csv
##     jobid: 13
##     reason: Missing output files: output/MR/EN-RAGE.csv
##     wildcards: protein=EN-RAGE
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 13:04:22 2022]
## rule MR_analysis:
##     input: input/MR.csv, input/ld/CX3CL1.ld, input/ld/CX3CL1.snplist
##     output: output/MR/CX3CL1.csv
##     jobid: 10
##     reason: Missing output files: output/MR/CX3CL1.csv
##     wildcards: protein=CX3CL1
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 13:04:22 2022]
## rule MR_analysis:
##     input: input/MR.csv, input/ld/SIRT2.ld, input/ld/SIRT2.snplist
##     output: output/MR/SIRT2.csv
##     jobid: 24
##     reason: Missing output files: output/MR/SIRT2.csv
##     wildcards: protein=SIRT2
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 13:04:22 2022]
## rule MR_analysis:
##     input: input/MR.csv, input/ld/IL-6.ld, input/ld/IL-6.snplist
##     output: output/MR/IL-6.csv
##     jobid: 17
##     reason: Missing output files: output/MR/IL-6.csv
##     wildcards: protein=IL-6
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 13:04:22 2022]
## rule MR_analysis:
##     input: input/MR.csv, input/ld/Beta-NGF.ld, input/ld/Beta-NGF.snplist
##     output: output/MR/Beta-NGF.csv
##     jobid: 3
##     reason: Missing output files: output/MR/Beta-NGF.csv
##     wildcards: protein=Beta-NGF
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 13:04:22 2022]
## rule MR_analysis:
##     input: input/MR.csv, input/ld/CCL4.ld, input/ld/CCL4.snplist
##     output: output/MR/CCL4.csv
##     jobid: 7
##     reason: Missing output files: output/MR/CCL4.csv
##     wildcards: protein=CCL4
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 13:04:22 2022]
## rule MR_analysis:
##     input: input/MR.csv, input/ld/MMP-10.ld, input/ld/MMP-10.snplist
##     output: output/MR/MMP-10.csv
##     jobid: 21
##     reason: Missing output files: output/MR/MMP-10.csv
##     wildcards: protein=MMP-10
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 13:04:22 2022]
## rule obs_analysis:
##     input: input/Obs.csv
##     output: output/Obs.csv
##     jobid: 1
##     reason: Missing output files: output/Obs.csv
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 13:04:22 2022]
## rule MR_analysis:
##     input: input/MR.csv, input/ld/FGF-23.ld, input/ld/FGF-23.snplist
##     output: output/MR/FGF-23.csv
##     jobid: 14
##     reason: Missing output files: output/MR/FGF-23.csv
##     wildcards: protein=FGF-23
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 13:04:22 2022]
## rule MR_analysis:
##     input: input/MR.csv, input/ld/VEGF-A.ld, input/ld/VEGF-A.snplist
##     output: output/MR/VEGF-A.csv
##     jobid: 28
##     reason: Missing output files: output/MR/VEGF-A.csv
##     wildcards: protein=VEGF-A
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 13:04:22 2022]
## rule MR_analysis:
##     input: input/MR.csv, input/ld/CD40.ld, input/ld/CD40.snplist
##     output: output/MR/CD40.csv
##     jobid: 8
##     reason: Missing output files: output/MR/CD40.csv
##     wildcards: protein=CD40
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 13:04:22 2022]
## rule MR_analysis:
##     input: input/MR.csv, input/ld/CXCL1.ld, input/ld/CXCL1.snplist
##     output: output/MR/CXCL1.csv
##     jobid: 11
##     reason: Missing output files: output/MR/CXCL1.csv
##     wildcards: protein=CXCL1
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 13:04:22 2022]
## rule MR_analysis:
##     input: input/MR.csv, input/ld/TNFSF14.ld, input/ld/TNFSF14.snplist
##     output: output/MR/TNFSF14.csv
##     jobid: 25
##     reason: Missing output files: output/MR/TNFSF14.csv
##     wildcards: protein=TNFSF14
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 13:04:22 2022]
## rule MR_analysis:
##     input: input/MR.csv, input/ld/CASP-8.ld, input/ld/CASP-8.snplist
##     output: output/MR/CASP-8.csv
##     jobid: 4
##     reason: Missing output files: output/MR/CASP-8.csv
##     wildcards: protein=CASP-8
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 13:04:22 2022]
## rule MR_analysis:
##     input: input/MR.csv, input/ld/IL-8.ld, input/ld/IL-8.snplist
##     output: output/MR/IL-8.csv
##     jobid: 18
##     reason: Missing output files: output/MR/IL-8.csv
##     wildcards: protein=IL-8
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 13:04:22 2022]
## rule MR_analysis:
##     input: input/MR.csv, input/ld/OPG.ld, input/ld/OPG.snplist
##     output: output/MR/OPG.csv
##     jobid: 22
##     reason: Missing output files: output/MR/OPG.csv
##     wildcards: protein=OPG
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 13:04:22 2022]
## rule MR_analysis:
##     input: input/MR.csv, input/ld/HGF.ld, input/ld/HGF.snplist
##     output: output/MR/HGF.csv
##     jobid: 15
##     reason: Missing output files: output/MR/HGF.csv
##     wildcards: protein=HGF
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 13:04:22 2022]
## rule aggregate_MR:
##     input: output/MR/Beta-NGF.csv, output/MR/CASP-8.csv, output/MR/CCL20.csv, output/MR/CCL3.csv, output/MR/CCL4.csv, output/MR/CD40.csv, output/MR/CSF-1.csv, output/MR/CX3CL1.csv, output/MR/CXCL1.csv, output/MR/CXCL6.csv, output/MR/EN-RAGE.csv, output/MR/FGF-23.csv, output/MR/HGF.csv, output/MR/IL-18.csv, output/MR/IL-6.csv, output/MR/IL-8.csv, output/MR/MCP-1.csv, output/MR/MMP-1.csv, output/MR/MMP-10.csv, output/MR/OPG.csv, output/MR/SCF.csv, output/MR/SIRT2.csv, output/MR/TNFSF14.csv, output/MR/TRAIL.csv, output/MR/TRANCE.csv, output/MR/VEGF-A.csv
##     output: output/MR.csv
##     jobid: 2
##     reason: Missing output files: output/MR.csv; Input files updated by another job: output/MR/CCL4.csv, output/MR/CCL20.csv, output/MR/FGF-23.csv, output/MR/TRAIL.csv, output/MR/CXCL1.csv, output/MR/CSF-1.csv, output/MR/IL-18.csv, output/MR/TNFSF14.csv, output/MR/IL-8.csv, output/MR/CD40.csv, output/MR/CX3CL1.csv, output/MR/EN-RAGE.csv, output/MR/Beta-NGF.csv, output/MR/SIRT2.csv, output/MR/OPG.csv, output/MR/TRANCE.csv, output/MR/CCL3.csv, output/MR/MMP-1.csv, output/MR/SCF.csv, output/MR/VEGF-A.csv, output/MR/CASP-8.csv, output/MR/CXCL6.csv, output/MR/IL-6.csv, output/MR/HGF.csv, output/MR/MMP-10.csv, output/MR/MCP-1.csv
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## 
## [Fri Jul 29 13:04:22 2022]
## localrule all:
##     input: output/Obs.csv, output/MR.csv
##     jobid: 0
##     reason: Input files updated by another job: output/MR.csv, output/Obs.csv
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
## [Fri Jul 29 13:04:23 2022]
## rule MR_analysis:
##     input: input/MR.csv, input/ld/CCL20.ld, input/ld/CCL20.snplist
##     output: output/MR/CCL20.csv
##     jobid: 5
##     reason: Missing output files: output/MR/CCL20.csv
##     wildcards: protein=CCL20
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 13:04:23 2022]
## rule MR_analysis:
##     input: input/MR.csv, input/ld/MCP-1.ld, input/ld/MCP-1.snplist
##     output: output/MR/MCP-1.csv
##     jobid: 19
##     reason: Missing output files: output/MR/MCP-1.csv
##     wildcards: protein=MCP-1
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 13:04:23 2022]
## rule MR_analysis:
##     input: input/MR.csv, input/ld/TRAIL.ld, input/ld/TRAIL.snplist
##     output: output/MR/TRAIL.csv
##     jobid: 26
##     reason: Missing output files: output/MR/TRAIL.csv
##     wildcards: protein=TRAIL
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 13:04:23 2022]
## rule MR_analysis:
##     input: input/MR.csv, input/ld/CXCL6.ld, input/ld/CXCL6.snplist
##     output: output/MR/CXCL6.csv
##     jobid: 12
##     reason: Missing output files: output/MR/CXCL6.csv
##     wildcards: protein=CXCL6
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 13:04:23 2022]
## rule MR_analysis:
##     input: input/MR.csv, input/ld/CCL3.ld, input/ld/CCL3.snplist
##     output: output/MR/CCL3.csv
##     jobid: 6
##     reason: Missing output files: output/MR/CCL3.csv
##     wildcards: protein=CCL3
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 13:04:23 2022]
## rule MR_analysis:
##     input: input/MR.csv, input/ld/MMP-1.ld, input/ld/MMP-1.snplist
##     output: output/MR/MMP-1.csv
##     jobid: 20
##     reason: Missing output files: output/MR/MMP-1.csv
##     wildcards: protein=MMP-1
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 13:04:23 2022]
## rule MR_analysis:
##     input: input/MR.csv, input/ld/TRANCE.ld, input/ld/TRANCE.snplist
##     output: output/MR/TRANCE.csv
##     jobid: 27
##     reason: Missing output files: output/MR/TRANCE.csv
##     wildcards: protein=TRANCE
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 13:04:23 2022]
## rule MR_analysis:
##     input: input/MR.csv, input/ld/EN-RAGE.ld, input/ld/EN-RAGE.snplist
##     output: output/MR/EN-RAGE.csv
##     jobid: 13
##     reason: Missing output files: output/MR/EN-RAGE.csv
##     wildcards: protein=EN-RAGE
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 13:04:23 2022]
## rule MR_analysis:
##     input: input/MR.csv, input/ld/CCL4.ld, input/ld/CCL4.snplist
##     output: output/MR/CCL4.csv
##     jobid: 7
##     reason: Missing output files: output/MR/CCL4.csv
##     wildcards: protein=CCL4
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 13:04:23 2022]
## rule MR_analysis:
##     input: input/MR.csv, input/ld/MMP-10.ld, input/ld/MMP-10.snplist
##     output: output/MR/MMP-10.csv
##     jobid: 21
##     reason: Missing output files: output/MR/MMP-10.csv
##     wildcards: protein=MMP-10
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 13:04:23 2022]
## rule obs_analysis:
##     input: input/Obs.csv
##     output: output/Obs.csv
##     jobid: 1
##     reason: Missing output files: output/Obs.csv
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 13:04:23 2022]
## rule MR_analysis:
##     input: input/MR.csv, input/ld/FGF-23.ld, input/ld/FGF-23.snplist
##     output: output/MR/FGF-23.csv
##     jobid: 14
##     reason: Missing output files: output/MR/FGF-23.csv
##     wildcards: protein=FGF-23
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 13:04:23 2022]
## rule MR_analysis:
##     input: input/MR.csv, input/ld/VEGF-A.ld, input/ld/VEGF-A.snplist
##     output: output/MR/VEGF-A.csv
##     jobid: 28
##     reason: Missing output files: output/MR/VEGF-A.csv
##     wildcards: protein=VEGF-A
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 13:04:23 2022]
## rule MR_analysis:
##     input: input/MR.csv, input/ld/OPG.ld, input/ld/OPG.snplist
##     output: output/MR/OPG.csv
##     jobid: 22
##     reason: Missing output files: output/MR/OPG.csv
##     wildcards: protein=OPG
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 13:04:23 2022]
## rule MR_analysis:
##     input: input/MR.csv, input/ld/CD40.ld, input/ld/CD40.snplist
##     output: output/MR/CD40.csv
##     jobid: 8
##     reason: Missing output files: output/MR/CD40.csv
##     wildcards: protein=CD40
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 13:04:23 2022]
## rule MR_analysis:
##     input: input/MR.csv, input/ld/HGF.ld, input/ld/HGF.snplist
##     output: output/MR/HGF.csv
##     jobid: 15
##     reason: Missing output files: output/MR/HGF.csv
##     wildcards: protein=HGF
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 13:04:23 2022]
## rule MR_analysis:
##     input: input/MR.csv, input/ld/CSF-1.ld, input/ld/CSF-1.snplist
##     output: output/MR/CSF-1.csv
##     jobid: 9
##     reason: Missing output files: output/MR/CSF-1.csv
##     wildcards: protein=CSF-1
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 13:04:23 2022]
## rule MR_analysis:
##     input: input/MR.csv, input/ld/SCF.ld, input/ld/SCF.snplist
##     output: output/MR/SCF.csv
##     jobid: 23
##     reason: Missing output files: output/MR/SCF.csv
##     wildcards: protein=SCF
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 13:04:23 2022]
## rule MR_analysis:
##     input: input/MR.csv, input/ld/IL-18.ld, input/ld/IL-18.snplist
##     output: output/MR/IL-18.csv
##     jobid: 16
##     reason: Missing output files: output/MR/IL-18.csv
##     wildcards: protein=IL-18
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 13:04:23 2022]
## rule MR_analysis:
##     input: input/MR.csv, input/ld/CX3CL1.ld, input/ld/CX3CL1.snplist
##     output: output/MR/CX3CL1.csv
##     jobid: 10
##     reason: Missing output files: output/MR/CX3CL1.csv
##     wildcards: protein=CX3CL1
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 13:04:23 2022]
## rule MR_analysis:
##     input: input/MR.csv, input/ld/SIRT2.ld, input/ld/SIRT2.snplist
##     output: output/MR/SIRT2.csv
##     jobid: 24
##     reason: Missing output files: output/MR/SIRT2.csv
##     wildcards: protein=SIRT2
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## 
## [Fri Jul 29 13:04:23 2022]
## rule MR_analysis:
##     input: input/MR.csv, input/ld/IL-6.ld, input/ld/IL-6.snplist
##     output: output/MR/IL-6.csv
##     jobid: 17
##     reason: Missing output files: output/MR/IL-6.csv
##     wildcards: protein=IL-6
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## 
## [Fri Jul 29 13:04:23 2022]
## rule MR_analysis:
##     input: input/MR.csv, input/ld/Beta-NGF.ld, input/ld/Beta-NGF.snplist
##     output: output/MR/Beta-NGF.csv
##     jobid: 3
##     reason: Missing output files: output/MR/Beta-NGF.csv
##     wildcards: protein=Beta-NGF
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## 
## [Fri Jul 29 13:04:23 2022]
## rule MR_analysis:
##     input: input/MR.csv, input/ld/CXCL1.ld, input/ld/CXCL1.snplist
##     output: output/MR/CXCL1.csv
##     jobid: 11
##     reason: Missing output files: output/MR/CXCL1.csv
##     wildcards: protein=CXCL1
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## 
## [Fri Jul 29 13:04:23 2022]
## rule MR_analysis:
##     input: input/MR.csv, input/ld/TNFSF14.ld, input/ld/TNFSF14.snplist
##     output: output/MR/TNFSF14.csv
##     jobid: 25
##     reason: Missing output files: output/MR/TNFSF14.csv
##     wildcards: protein=TNFSF14
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## 
## [Fri Jul 29 13:04:23 2022]
## rule MR_analysis:
##     input: input/MR.csv, input/ld/CASP-8.ld, input/ld/CASP-8.snplist
##     output: output/MR/CASP-8.csv
##     jobid: 4
##     reason: Missing output files: output/MR/CASP-8.csv
##     wildcards: protein=CASP-8
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## 
## [Fri Jul 29 13:04:23 2022]
## rule MR_analysis:
##     input: input/MR.csv, input/ld/IL-8.ld, input/ld/IL-8.snplist
##     output: output/MR/IL-8.csv
##     jobid: 18
##     reason: Missing output files: output/MR/IL-8.csv
##     wildcards: protein=IL-8
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.2 ──
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
## ✔ ggplot2 3.3.6     ✔ purrr   0.3.4
## ✔ tibble  3.1.8     ✔ dplyr   1.0.9
## ✔ tidyr   1.2.0     ✔ stringr 1.4.0
## ✔ readr   2.1.2     ✔ forcats 0.5.1
## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.2 ──
## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.2 ──
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
## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
## ✖ dplyr::between()   masks data.table::between()
## ✖ dplyr::filter()    masks stats::filter()
## ✖ dplyr::first()     masks data.table::first()
## ✖ dplyr::lag()       masks stats::lag()
## ✖ dplyr::last()      masks data.table::last()
## ✖ purrr::transpose() masks data.table::transpose()
## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
## ✖ dplyr::filter() masks stats::filter()
## ✖ dplyr::lag()    masks stats::lag()
## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.2 ──
## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.2 ──
## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.2 ──
## Loading required package: Matrix
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
## Rows: 124 Columns: 5
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
## ℹ The warning occurred in group 26: Protein = "SIRT2". 
## 3: `cols` is now required when using unnest().
## Please use `cols = c(meta)` 
## [Fri Jul 29 13:04:28 2022]
## Finished job 1.
## 1 of 29 steps (3%) done
## Method requires data on >2 variants.Method requires data on >2 variants.Method requires data on >2 variants.Method requires data on >2 variants.Method requires data on >2 variants.Method requires data on >2 variants.Method requires data on >2 variants.Method requires data on >2 variants.Method requires data on >2 variants.Method requires data on >2 variants.Method requires data on >2 variants.Method requires data on >2 variants.Method requires data on >2 variants.Method requires data on >2 variants.Method requires data on >2 variants.Method requires data on >2 variants.Method requires data on >2 variants.Method requires data on >2 variants.Method requires data on >2 variants.Method requires data on >2 variants.Method requires data on >2 variants.Method requires data on >2 variants.Method requires data on >2 variants.Method requires data on >2 variants.Method requires data on >2 variants.Method requires data on >2 variants.Method requires data on >2 variants.Method requires data on >2 variants.Method requires data on >2 variants.Method requires data on >2 variants.Method requires data on >2 variants.Method requires data on >2 variants.[Fri Jul 29 13:04:29 2022]
## Finished job 23.
## 2 of 29 steps (7%) done
## Method requires data on >2 variants.Method requires data on >2 variants.[Fri Jul 29 13:04:29 2022]
## Finished job 18.
## 3 of 29 steps (10%) done
## [Fri Jul 29 13:04:29 2022]
## Finished job 4.
## 4 of 29 steps (14%) done
## [Fri Jul 29 13:04:30 2022]
## Finished job 16.
## 5 of 29 steps (17%) done
## [Fri Jul 29 13:04:30 2022]
## Finished job 7.
## 6 of 29 steps (21%) done
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
## [Fri Jul 29 13:04:30 2022]
## Finished job 11.
## 7 of 29 steps (24%) done
## [Fri Jul 29 13:04:30 2022]
## Finished job 17.
## 8 of 29 steps (28%) done
## [Fri Jul 29 13:04:30 2022]
## Finished job 19.
## 9 of 29 steps (31%) done
## Warning messages:
## 1: Problem while computing `MR = map(data, run_MR_all, ldrho)`.
## ℹ NaNs produced
## ℹ The warning occurred in group 23: r2_thresh = 0.4, P_thresh = 0.01, Protein =
##   "TRANCE". 
## 2: Problem while computing `MR = map(data, run_MR_all, ldrho)`.
## ℹ NaNs produced
## ℹ The warning occurred in group 23: r2_thresh = 0.4, P_thresh = 0.01, Protein =
##   "TRANCE". 
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
## [Fri Jul 29 13:04:30 2022]
## Finished job 12.
## 10 of 29 steps (34%) done
## Warning messages:
## 1: Problem while computing `MR = map(data, run_MR_all, ldrho)`.
## ℹ NaNs produced
## ℹ The warning occurred in group 30: r2_thresh = 0.6, P_thresh = 1, Protein =
##   "CCL3". 
## 2: Problem while computing `MR = map(data, run_MR_all, ldrho)`.
## ℹ NaNs produced
## ℹ The warning occurred in group 30: r2_thresh = 0.6, P_thresh = 1, Protein =
##   "CCL3". 
## [Fri Jul 29 13:04:30 2022]
## Finished job 27.
## 11 of 29 steps (38%) done
## [Fri Jul 29 13:04:30 2022]
## Finished job 13.
## 12 of 29 steps (41%) done
## [Fri Jul 29 13:04:30 2022]
## Finished job 6.
## 13 of 29 steps (45%) done
## Warning message:
## Problem while computing `MR = map(data, run_MR_all, ldrho)`.
## ℹ NaNs produced
## ℹ The warning occurred in group 30: r2_thresh = 0.6, P_thresh = 1, Protein =
##   "HGF". 
## Method requires data on >2 variants.Warning messages:
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
## [Fri Jul 29 13:04:30 2022]
## Finished job 15.
## 14 of 29 steps (48%) done
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
## [Fri Jul 29 13:04:30 2022]
## Finished job 24.
## 15 of 29 steps (52%) done
## [Fri Jul 29 13:04:30 2022]
## Finished job 22.
## 16 of 29 steps (55%) done
## Warning message:
## Problem while computing `MR = map(data, run_MR_all, ldrho)`.
## ℹ NaNs produced
## ℹ The warning occurred in group 25: r2_thresh = 0.6, P_thresh = 5e-08, Protein
##   = "CX3CL1". 
## [Fri Jul 29 13:04:31 2022]
## Finished job 10.
## 17 of 29 steps (59%) done
## [Fri Jul 29 13:04:31 2022]
## Finished job 9.
## 18 of 29 steps (62%) done
## [Fri Jul 29 13:04:31 2022]
## Finished job 3.
## 19 of 29 steps (66%) done
## [Fri Jul 29 13:04:31 2022]
## Finished job 8.
## 20 of 29 steps (69%) done
## [Fri Jul 29 13:04:32 2022]
## Finished job 28.
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
## [Fri Jul 29 13:04:32 2022]
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
## [Fri Jul 29 13:04:32 2022]
## Finished job 5.
## 23 of 29 steps (79%) done
## [Fri Jul 29 13:04:33 2022]
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
## Warning message:
## Problem while computing `MR = map(data, run_MR_all, ldrho)`.
## ℹ NaNs produced
## ℹ The warning occurred in group 24: r2_thresh = 0.4, P_thresh = 1, Protein =
##   "MMP-10". 
## [Fri Jul 29 13:04:33 2022]
## Finished job 26.
## 25 of 29 steps (86%) done
## [Fri Jul 29 13:04:33 2022]
## Finished job 21.
## 26 of 29 steps (90%) done
## [Fri Jul 29 13:04:35 2022]
## Finished job 25.
## 27 of 29 steps (93%) done
## Select jobs to execute...
## 
## [Fri Jul 29 13:04:35 2022]
## rule aggregate_MR:
##     input: output/MR/Beta-NGF.csv, output/MR/CASP-8.csv, output/MR/CCL20.csv, output/MR/CCL3.csv, output/MR/CCL4.csv, output/MR/CD40.csv, output/MR/CSF-1.csv, output/MR/CX3CL1.csv, output/MR/CXCL1.csv, output/MR/CXCL6.csv, output/MR/EN-RAGE.csv, output/MR/FGF-23.csv, output/MR/HGF.csv, output/MR/IL-18.csv, output/MR/IL-6.csv, output/MR/IL-8.csv, output/MR/MCP-1.csv, output/MR/MMP-1.csv, output/MR/MMP-10.csv, output/MR/OPG.csv, output/MR/SCF.csv, output/MR/SIRT2.csv, output/MR/TNFSF14.csv, output/MR/TRAIL.csv, output/MR/TRANCE.csv, output/MR/VEGF-A.csv
##     output: output/MR.csv
##     jobid: 2
##     reason: Missing output files: output/MR.csv; Input files updated by another job: output/MR/MCP-1.csv, output/MR/VEGF-A.csv, output/MR/CCL3.csv, output/MR/MMP-10.csv, output/MR/EN-RAGE.csv, output/MR/CASP-8.csv, output/MR/CXCL1.csv, output/MR/CXCL6.csv, output/MR/HGF.csv, output/MR/TRAIL.csv, output/MR/TRANCE.csv, output/MR/CD40.csv, output/MR/IL-18.csv, output/MR/CX3CL1.csv, output/MR/SCF.csv, output/MR/OPG.csv, output/MR/Beta-NGF.csv, output/MR/FGF-23.csv, output/MR/CCL20.csv, output/MR/TNFSF14.csv, output/MR/CCL4.csv, output/MR/IL-6.csv, output/MR/CSF-1.csv, output/MR/SIRT2.csv, output/MR/MMP-1.csv, output/MR/IL-8.csv
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 13:04:35 2022]
## Finished job 2.
## 28 of 29 steps (97%) done
## Select jobs to execute...
## 
## [Fri Jul 29 13:04:35 2022]
## localrule all:
##     input: output/Obs.csv, output/MR.csv
##     jobid: 0
##     reason: Input files updated by another job: output/Obs.csv, output/MR.csv
##     resources: tmpdir=/rds/user/jhz22/hpc-work/work
## 
## [Fri Jul 29 13:04:35 2022]
## Finished job 0.
## 29 of 29 steps (100%) done
## Complete log: .snakemake/log/2022-07-29T130422.918957.snakemake.log
```

which gives `output`/`MR.csv` (MR results) and `Obs.csv` (observational results).
