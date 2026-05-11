#!/usr/bin/env bash
#SBATCH --job-name=pQTLtools
#SBATCH --account PETERS-SL3-CPU
#SBATCH --partition icelake-himem
#SBATCH --mem=28800
#SBATCH --time=12:00:00
#SBATCH --output=/home/jhz22/R/work/pQTLtools_%A.out
#SBATCH --error=/home/jhz22/R/work/pQTLtools_%A.err
#SBATCH --export=ALL

set -euo pipefail

. /etc/profile.d/modules.sh
module purge
module load rhel8/default-icl ceuadmin/R

export TMPDIR=/rds/user/jhz22/hpc-work/work
export IN_PKGDOWN=true
export R_PKGDOWN_BUILD_LLM=false

cd ~/pQTLtools

Rscript -e '
  suppressMessages({library(pkgdown);library(roxygen2);library(knitr)})
  knitr::knit("README.Rmd")
  roxygen2::roxygenise()
  pkgdown::build_site()
'

git add .
git commit -m "pkgdown rebuild $(date +%F)" || true
git push || true

du -hs --exclude .git --exclude docs
