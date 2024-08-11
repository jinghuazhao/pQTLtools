#!/usr/bin/bash

if [ "$(uname -n | sed 's/-[0-9]*$//')" == "login-q" ]; then
   module load ceuadmin/libssh/0.10.6-icelake
   module load ceuadmin/openssh/9.7p1-icelake
   module load bcftools/1.13/gcc/intel-oneapi-mkl/vi7yswnv
fi
module load ceuadmin/R
module load mono/5.0.1.1
module load texlive

export version=$(awk '/Version/{print $2}' ~/pQTLtools/DESCRIPTION)
R CMD build --resave-data --compact-vignettes=both pQTLtools
R CMD INSTALL pQTLtools_${version}.tar.gz

R CMD check --as-cran pQTLtools_${version}.tar.gz
rm -rf ${HOME}/R/pQTLtools.Rcheck
mv pQTLtools_${version}.tar.gz pQTLtools.Rcheck R
