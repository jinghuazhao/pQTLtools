# UniProt IDs to others

A function which converts UniProt IDs to others.

## Usage

``` r
uniprot2ids(uniprotid = "ACC+ID", to, query)
```

## Arguments

-   uniprotid:

    Source IDs.

-   to:

    To IDs.

-   query:

    A query.

## Value

A UniProt-ID mapping

## Details

This function is based on the Python3 script from UniProt. See
https://www.uniprot.org/help/api_idmapping

## Note

Adapted from script by UniProt

## Examples

``` r
if (FALSE) { # \dontrun{
uniprotid <- "ACC+ID"
to <- "CHEMBL_ID"
query <- noquote(gap.datasets::inf1[["uniprot"]])
query <- paste(query,collapse=" ")
r <- pQTLtools::uniprot2ids(uniprotid,to,query)
cat(r,file="INF1.merge.chembl")
} # }
```
