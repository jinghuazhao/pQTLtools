#!/usr/bin/bash

function setup()
{
Rscript -e '
  suppressMessages(library(pkgdown))
  # DESCRIPTION
  {
    knitr::knit("README.Rmd")
    roxygen2::roxygenise()
    pkgdown::build_site()
  }
  pkgdown::build_site_github_pages()
# usethis::use_github_action("pkgdown", save_as = "R-CMD-check.yaml", ref = NULL, ignore = TRUE, open = FALSE)
# clean_site(); init_site(); build_home(); build_news(); build_articles(); build_reference(); build_search()
'
# devtools::build_rmd() is equivalent but limited, so knitr::knit/pandoc are better options.
# pandoc README.md --citeproc --mathjax -s --self-contained -o index.html
}

cp vignettes/ACE.js vignettes/lz.htm docs/articles
if [ -d docs/articles/data ]; then rm -rf docs/articles/data; cp -r vignettes/data docs/articles; fi
for d in bioconductor esse pQTLtools spectrum
do
    if [ -d vignettes/${d} ]; then
       rm -rf docs/articles/${d}
       mv vignettes/${d} docs/articles/
    fi
done
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
