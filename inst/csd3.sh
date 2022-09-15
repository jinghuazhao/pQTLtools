#!/usr/bin/bash

module load gcc/6 texlive

Rscript -e 'setwd("~/pQTLtools");devtools::document()'

export version=0.2
R CMD build --resave-data --compact-vignettes=both pQTLtools
R CMD INSTALL pQTLtools_${version}.tar.gz

R CMD check --as-cran pQTLtools_${version}.tar.gz
rm -rf ${HOME}/R/pQTLtools.Rcheck
mv pQTLtools_${version}.tar.gz pQTLtools.Rcheck R
