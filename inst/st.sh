#!/usr/bin/bash

function build()
{
  R --no-save -q <<\ \ END
    library(pkgdown)
  # keep as appropriate
    clean_site()
    build_home()
    build_news()
    build_articles()
    build_reference()
    build_site()
  END

  mv vignettes/gap docs/articles/
  mv vignettes/pQTLtools docs/articles/
}

for f in .github .gitignore .Rbuildignore .Rinstignore .travis.yml \
         data/ DESCRIPTION docs/ INDEX inst/ LICENSE LICENSE.md man/ NAMESPACE NEWS.md pkgdown/ R/ README.md vignettes/
do
  echo adding ${f}
  git add -f ${f}
  git commit -m "${f}"
done

git push
