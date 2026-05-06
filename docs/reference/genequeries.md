# phenoscanner genequeries in batches

R/phenoscanner only allows for certain number of items supplied. This
simple function return a large number of calls in batches as well as
generating SNPIDs.

## Usage

``` r
genequeries(
  genelist,
  catalogue = "pQTL",
  proxies = "EUR",
  p = 5e-08,
  r2 = 0.8,
  build = 37,
  wait = TRUE
)
```

## Arguments

-   genelist:

    a list of SNPs.

-   catalogue:

    "None","eQTL","mQTL","methQTL","pQTL","GWAS".

-   proxies:

    "None", "AFR","AMR","EAS","EUR","SAS".

-   p:

    p value threshold.

-   r2:

    r2 for LD.

-   build:

    37, 38.

-   wait:

    a flag to wait for 1hr for every 50 genes.

## Value

The returned value is a list containing genes and results.

## Details

Batches are generated and queries are combined into one.

## See also

`phenoscanner`

## Examples

``` r
if (FALSE) { # \dontrun{
# single gene
  genequeries("TNFRSF11B")
} # }
```
