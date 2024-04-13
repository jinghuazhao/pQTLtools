#!/usr/bin/bash

function build()
{
  Rscript -e 'devtools::build_readme();devtools::document()'
  Rscript -e '
    usethis::use_github_action("pkgdown", save_as = "R-CMD-check.yaml", ref = NULL, ignore = TRUE, open = FALSE)
    libs <- c("gwasvcf","rtracklayer","VariantAnnotation")
    invisible(suppressMessages(lapply(libs, require, character.only = TRUE)))
    library(pkgdown)
#   clean_site(); init_site(); build_home(); build_news(); build_articles(); build_reference(); build_search()
    build_site()
  '
# pandoc README.md --citeproc --mathjax -s --self-contained -o index.html
}

build

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

for f in .github .gitignore .Rbuildignore .Rinstignore .travis.yml \
         DESCRIPTION INDEX LICENSE LICENSE.md NAMESPACE NEWS.md R/ README.* \
         docs/ inst/ man/ pkgdown/ vignettes/
do
  echo adding ${f}
  git add ${f}
  git commit -m "${f}"
done

git push
du -hs --exclude .git --exclude docs
