#!/usr/bin/bash

module load ceuadmin/R
module load mono/5.0.1.1
module load texlive

export version=$(awk '/Version/{print $2}' ~/pQTLtools/DESCRIPTION)
R CMD build --resave-data --compact-vignettes=both pQTLtools
R CMD INSTALL pQTLtools_${version}.tar.gz

R CMD check --as-cran pQTLtools_${version}.tar.gz
rm -rf ${HOME}/R/pQTLtools.Rcheck
mv pQTLtools_${version}.tar.gz pQTLtools.Rcheck R
