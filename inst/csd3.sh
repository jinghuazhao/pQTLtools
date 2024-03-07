#!/usr/bin/bash

module load gcc/6 texlive ceuadmin/R

export version=$(awk '/Version/{print $2}' ~/pQTLtools/DESCRIPTION)
R CMD build --resave-data --compact-vignettes=both pQTLtools
R CMD INSTALL pQTLtools_${version}.tar.gz

R CMD check --as-cran pQTLtools_${version}.tar.gz
rm -rf ${HOME}/R/pQTLtools.Rcheck
mv pQTLtools_${version}.tar.gz pQTLtools.Rcheck R
