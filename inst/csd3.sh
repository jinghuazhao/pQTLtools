#!/usr/bin/bash

module load gcc/6
module load pcre/8.38
module load texlive

R CMD build --resave-data --compact-vignettes pQTLtools
R CMD check --as-cran pQTLtools_0.1.tar.gz
R CMD INSTALL pQTLtools_0.1.tar.gz

# also work:
# R CMD INSTALL pQTLtools

rm -rf pQTLtools/pQTLtools.Rcheck
mv pQTLtools_0.1.tar.gz R
