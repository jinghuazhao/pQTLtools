# Limit of detection analysis for ExpressionSet objects

Computes the percentage of values below a feature-specific lower limit
of detection (LOD) in an `ExpressionSet`. Commonly used in proteomic
quality assessment.

## Usage

``` r
get.prop.below.LLOD(eset, flagged = c("OUT", "IN"))
```

## Arguments

-   eset:

    An `ExpressionSet` object containing expression values and feature
    metadata. Must include a numeric column `lod.max` in `fData(eset)`.

-   flagged:

    Character string indicating whether flagged samples should be
    removed. One of `"OUT"` (remove flagged samples) or `"IN"` (retain
    all samples).

## Value

An `ExpressionSet` with an added feature annotation column:

-   pc.belowLOD.new:

    Percentage of values below the LOD per feature.

## Details

The function compares expression values from `exprs(eset)` to
feature-specific LOD thresholds stored in `fData(eset)$lod.max`.
Computation is vectorised using `sweep()` for efficiency and consistency
with Bioconductor conventions.

Flagged samples (if present in `pData(eset)$Flagged`) can optionally be
removed.

## Examples

``` r
if (FALSE) { # \dontrun{
suppressMessages(library(Biobase))
data(sample.ExpressionSet, package = "Biobase")
exampleSet <- sample.ExpressionSet
set.seed(1)
probs <- runif(nrow(exampleSet))
Biobase::fData(exampleSet)$lod.max <- vapply(seq_len(nrow(exampleSet)), function(i)
           quantile(Biobase::exprs(exampleSet)[i, ], probs = probs[i], na.rm = TRUE),
           numeric(1))
lod <- get.prop.below.LLOD(exampleSet)
fd <- Biobase::fData(lod)
fd <- fd[order(fd$pc.belowLOD.new, decreasing = TRUE), ]
knitr::kable(head(fd))
plot(fd$pc.belowLOD.new, main = "Random quantile cut off", ylab = "% below LOD")
} # }
```
