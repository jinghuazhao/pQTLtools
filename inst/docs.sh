#!/usr/bin/bash

Rscript -e '
  devtools::build_readme();devtools::document()
  pkgdown::build_site()
# usethis::use_github_action("pkgdown", save_as = "R-CMD-check.yaml", ref = NULL, ignore = TRUE, open = FALSE)
# clean_site(); init_site(); build_home(); build_news(); build_articles(); build_reference(); build_search()
'
# pandoc README.md --citeproc --mathjax -s --self-contained -o index.html

if [ -d vignettes/bioconductor ]; then
   rm -rf docs/articles/bioconductor
   mv vignettes/bioconductor docs/articles/
fi
if [ -d vignettes/es ]; then
   rm -rf docs/articles/es
   mv vignettes/es docs/articles/
fi
if [ -d vignettes/pQTLtools ]; then
   rm -rf docs/articles/pQTLtools
   mv vignettes/pQTLtools docs/articles/
fi
if [ -f vignettes/fig2d.html ]; then
   mv vignettes/fig2d.html docs/articles/
fi

if [ -f vignettes/fig3d.html ]; then
   mv vignettes/fig3d.html docs/articles/
fi

if [ "$(uname -n | sed 's/-[0-9]*$//')" == "login-q" ]; then
   module load ceuadmin/libssh/0.10.6-icelake
   module load ceuadmin/openssh/9.7p1-icelake
fi

for f in .gitignore .Rbuildignore .Rinstignore .travis.yml \
         DESCRIPTION LICENSE LICENSE.md NAMESPACE NEWS.md R/ README.Rmd README.md \
         docs/ inst/ man/ pkgdown/ vignettes/
do
  echo adding ${f}
  git add ${f}
  git commit -m  --no-verify "${f}"
done

git push
du -hs --exclude .git --exclude docs
