#!/usr/bin/bash

function setup()
{
Rscript -e '
  devtools::build_readme()
  roxygen2::roxygenise()
  pkgdown::build_site_github_pages()
# usethis::use_github_action("pkgdown", save_as = "R-CMD-check.yaml", ref = NULL, ignore = TRUE, open = FALSE)
# clean_site(); init_site(); build_home(); build_news(); build_articles(); build_reference(); build_search()
'
# pandoc README.md --citeproc --mathjax -s --self-contained -o index.html
}

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
  git commit -m  --no-verify "${f}"
done

git push
du -hs --exclude .git --exclude docs
