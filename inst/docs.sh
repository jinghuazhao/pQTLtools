#!/usr/bin/bash

function build()
{
  Rscript -e "devtools::build_readme()"
  Rscript -e "
    libs <- c('gwasvcf','rtracklayer','VariantAnnotation')
    invisible(suppressMessages(lapply(libs, require, character.only = TRUE)))
    devtools::document()
    library(pkgdown)
    clean_site()
    init_site()
    build_home()
    build_news()
    build_articles()
    build_reference()
    build_search()
    build_site()
  "
}

build

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
         data/ DESCRIPTION docs/ INDEX inst/ LICENSE LICENSE.md man/ NAMESPACE NEWS.md pkgdown/ R/ README.* vignettes/
do
  echo adding ${f}
  git add -f ${f}
  git commit -m "${f}"
done

git push -f
