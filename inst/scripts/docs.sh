#!/usr/bin/env bash

#SBATCH --job-name=_pQTLtools
#SBATCH --account PETERS-SL3-CPU
#SBATCH --partition icelake-himem
#SBATCH --mem=28800
#SBATCH --time=12:00:00
#SBATCH --error=/home/jhz22/R/work/_pQTLtools_%A_%a.err
#SBATCH --output=/home/jhz22/R/work/_pQTLtools_%A_%a.out
#SBATCH --export ALL

. /etc/profile.d/modules.sh
module purge
module load rhel8/default-icl
module load ceuadmin/R

export TMPDIR=/rds/user/jhz22/hpc-work/work

cd ~/pQTLtools
Rscript -e '
  suppressMessages(library(pkgdown))
  {
    knitr::knit("README.Rmd")
    roxygen2::roxygenise()
    pkgdown::build_site()
  }
# devtools::build_rmd() is equivalent but limited compared to knitr::knit/pandoc.
# pkgdown::build_site_github_pages()
# usethis::use_github_action("pkgdown", save_as = "R-CMD-check.yaml", ref = NULL, ignore = TRUE, open = FALSE)
# clean_site(); init_site(); build_home(); build_news(); build_articles(); build_reference(); build_search()
'
# pandoc README.md --citeproc --mathjax -s --self-contained -o index.html

cp -r vignettes/data vignettes/lz.html vignettes/stack.html docs/articles
for d in bioconductor pQTLtools snakemake
do
    if [ -d vignettes/${d} ]; then
       rm -rf docs/articles/${d}
       mv vignettes/${d} docs/articles/
    fi
done
if [ -f vignettes/ebi-a-GCST007432.vcf.gz.tbi ]; then rm vignettes/ebi-a-GCST007432.vcf.gz.tbi; fi
for f in fig2d fig3d heatmaply mcpca3d
do
    if [ -f vignettes/${f}.html ]; then mv vignettes/${f}.html docs/articles/; fi
done

for f in .gitignore .Rbuildignore .Rinstignore .travis.yml \
         DESCRIPTION LICENSE LICENSE.md NAMESPACE NEWS.md README.Rmd README.md \
         R/ data/ docs/ inst/ man/ pkgdown/ vignettes/
do
  echo adding ${f}
  git add ${f}
  git commit --no-verify -m "${f}"
done

git push
du -hs --exclude .git --exclude docs
