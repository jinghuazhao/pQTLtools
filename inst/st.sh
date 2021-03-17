#!/usr/bin/bash

Rscript -e "library(pkgdown);build_articles();build_news();build_reference()"

for f in .github .gitignore .Rbuildignore .Rinstignore .travis.yml \
         data/ DESCRIPTION docs/ INDEX inst/ LICENSE LICENSE.md man/ NAMESPACE NEWS.md pkgdown/ R/ README.md vignettes/
do
  echo adding ${f}
  git add -f ${f}
  git commit -m "${f}"
done

git push
