# Locus novelty check using LD overlap

This function assesses whether query loci are novel or match known loci
by combining genomic proximity (± flanking window) and linkage
disequilibrium (LD).

## Usage

``` r
novelty_check(
  known_loci,
  query_loci,
  ldops = NULL,
  flanking = 1e+06,
  pop = "EUR",
  verbose = TRUE
)
```

## Arguments

-   known_loci:

    Data.frame of known/published loci. Must contain columns: chr, pos,
    uniprot, rsid, prot.

-   query_loci:

    Data.frame of query loci to evaluate. Must contain columns: chr,
    pos, uniprot, rsid, prot.

-   ldops:

    Optional list specifying local LD computation:

    bfile

    :   PLINK binary prefix (bed/bim/fam)

    plink

    :   Path to PLINK executable

-   flanking:

    Genomic window (bp) around query loci used for overlap. Default is
    1e6 (±1 Mb).

-   pop:

    1000 Genomes population code (e.g., "EUR") used when ldops = NULL.

-   verbose:

    Logical; if TRUE prints missing LD variants.

## Value

A data.frame with paired known/query loci and LD r^2 values.

## Details

It supports:

-   1000 Genomes LD reference via ieugwasr

-   Local PLINK reference panels via ld_matrix_local()

The function:

1.  Finds overlapping loci within a flanking distance

2.  Matches loci by gene/protein (uniprot)

3.  Computes LD (r and r^2) between known and query variants

4.  Returns per-pair LD values for downstream novelty/replication
    assessment

## Examples

``` r
if (FALSE) { # \dontrun{
# 1000G mode
novelty_check(known_loci, query_loci)

# Local PLINK mode
novelty_check(
  known_loci,
  query_loci,
  ldops = list(
    bfile = "/path/interval.imputed.olink.chr_3",
    plink = "/path/plink"
  )
)
} # }
```
