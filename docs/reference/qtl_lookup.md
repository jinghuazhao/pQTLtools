# QTL lookup with LD-based proxy selection

This function performs QTL lookup by integrating MR results with GWAS
summary statistics and identifying proxy SNPs based on linkage
disequilibrium (LD).

## Usage

``` r
qtl_lookup(
  d,
  dat,
  panel = c("1000Genomes", "local"),
  p_threshold = 0.001,
  r2_threshold = 0.8,
  pop = "EUR",
  plink_bin = NULL,
  max_snps = 500,
  xlsx = NULL,
  verbose = TRUE
)
```

## Arguments

-   d:

    Directory containing GWAS files.

-   dat:

    Data frame of MR results. Must contain: protein, id, pqtl, qtl,
    p_qtl, file_gwas.

-   panel:

    LD reference type: "1000Genomes" or "local".

-   p_threshold:

    P-value cutoff for selecting candidate SNPs.

-   r2_threshold:

    LD r^2 threshold for locus classification.

-   pop:

    1000 Genomes population (default "EUR").

-   plink_bin:

    Path to PLINK executable (required for local mode).

-   max_snps:

    Maximum SNPs used in LD computation.

-   xlsx:

    Optional Excel output path.

-   verbose:

    Logical; show progress messages.

## Value

Updated dat with columns: proxy, p_proxy, rsq, classification.

## Details

It supports:

-   1000 Genomes LD reference (ieugwasr API)

-   Local PLINK reference panels (recommended for large-scale runs)

Loci are classified as:

-   same_locus (r^2 \>= threshold)

-   independent_locus (r^2 \< threshold)
