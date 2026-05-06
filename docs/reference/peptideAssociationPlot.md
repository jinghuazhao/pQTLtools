# peptide association plot

peptide association plot

## Usage

``` r
peptideAssociationPlot(protein, cistrans, chrlen = gap::hg19, disp = 85)
```

## Arguments

-   protein:

    protein name.

-   cistrans:

    peptide signals and classification.

-   chrlen:

    chromosome length (hg18, hg19, hg38).

-   disp:

    distance from x-axis.

## Value

None.

## Examples

``` r
if (FALSE) { # \dontrun{
  par(mfrow=c(2,1))
  protein <- "PROC"
  suffix <- "_dr"
  input <- paste0("~/Caprion/analysis/METAL",suffix,"/gz/",protein,suffix,".txt.gz")
  annotation <- paste0("~/Caprion/analysis/METAL",suffix,"/vep/",protein,suffix,".txt")
  reference <- file.path(find.package("pQTLtools"),"turboman",
                                      "turboman_hg19_reference_data.rda")
  pvalue_sign <- 5e-8
  plot_title <- protein
  pQTLtools::turboman(input, annotation, reference, pvalue_sign, plot_title)
  cistrans <- read.csv(paste0("~/pQTLtools/tests","/",protein,".cis.vs.trans"))
  load("~/pQTLtools/tests/PROC.rda")
  peptideAssociationPlot(protein,cistrans)
} # }
```
