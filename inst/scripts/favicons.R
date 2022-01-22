
# 1. generate logo.svg

library(gap)
svg("logo.svg")
INF <- Sys.getenv("INF")
d <- read.csv(file.path(INF,"work","INF1.merge.cis.vs.trans"),as.is=TRUE)
r <- pqtl2dplot(d)
dev.off()

# 2. build favicons

library(pkgdown)
build_favicons(overwrite=TRUE)
