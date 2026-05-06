# peptide-to-protein mapping

peptide-to-protein mapping

## Usage

``` r
peptideMapping(protein, batch = "ZWK", mm = 5)
```

## Arguments

-   protein:

    protein name.

-   batch:

    batch No. in the Caprion experiment.

-   mm:

    maximum number of mismatches.

## Value

a list containing mapping information.

## Examples

``` r
if (FALSE) { # \dontrun{
  batch <- "ZWK"
  load(paste0("~/Caprion/pilot/",batch,".rda"))
  PROC <- peptideMapping("PROC",mm=0)
} # }
```
