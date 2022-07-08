#!/usr/bin/bash

module load gcc/6
module load texlive

Rscript -e 'setwd("~/pQTLtools");devtools::document()'

R CMD build --resave-data --compact-vignettes=both pQTLtools
R CMD INSTALL pQTLtools_0.1.tar.gz

R CMD check --as-cran pQTLtools_0.1.tar.gz
rm -rf ${HOME}/R/pQTLtools.Rcheck
mv pQTLtools_0.1.tar.gz pQTLtools.Rcheck R
