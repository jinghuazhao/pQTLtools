#!/usr/bin/bash

for f in .gitignore .Rbuildignore .Rinstignore .travis.yml \
         data/ DESCRIPTION docs/ INDEX index.md inst/ LICENSE LICENSE.md man/ NAMESPACE NEWS.md pkgdown/ R/ st.sh vignettes/
do
  git add ${f}
  git commit -m "${f}"
done
