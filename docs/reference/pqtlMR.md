# Basic pQTL-MR analysis

This function takes exposure and outcome data as produced by
`format_data()` and used to perform MR analysis against a list of
outcomes; in the latter case it can be data from MR-Base, e.g.,
`outcome <- extract_outcome_data(snps=with(exposure,SNP),outcomes=c("ieu-a-7","ebi-a-GCST007432"))`.

## Usage

``` r
pqtlMR(
  exposure,
  outcome,
  mr_plot = FALSE,
  prefix = "pQTL-combined-",
  reverse = FALSE
)
```

## Arguments

-   exposure:

    exposure data.

-   outcome:

    the counterpart for outcome.

-   mr_plot:

    to produce plots.

-   prefix:

    a prefix for output files.

-   reverse:

    if TRUE, perform reverse MR.

## Value

No value is returned but several files.

## Details

Adapted from Zheng et al. (2020) , this function is analogous to
`run_TwoSampleMR()`.

## Note

Adapted from script by Jie Zheng.

## References

Zheng J, Haberland V, Baird D, Walker V, Haycock PC, Hurle MR,
Gutteridge A, Erola P, Liu Y, Luo S, Robinson J, Richardson TG, Staley
JR, Elsworth B, Burgess S, Sun BB, Danesh J, Runz H, Maranville JC,
Martin HM, Yarmolinsky J, Laurin C, Holmes MV, Liu JZ, Estrada K, Santos
R, McCarthy L, Waterworth D, Nelson MR, Smith GD, Butterworth AS, Hemani
G, Scott RA, Gaunt TR (2020). “Phenome-wide Mendelian randomization
mapping the influence of the plasma proteome on complex diseases.”
*Nature Genetics*, **52**(10), 1122-1131.
[doi:10.1038/s41588-020-0682-6](https://doi.org/10.1038/s41588-020-0682-6)
.

## Examples

``` r
fi <- file.path(find.package("pQTLtools",lib.loc=.libPaths()),"tests","Ins.csv")
exposure <- TwoSampleMR::format_data(read.csv(fi),type="exposure")
fo <- file.path(find.package("pQTLtools",lib.loc=.libPaths()),"tests","Out.csv")
outcome <- TwoSampleMR::format_data(read.csv(fo),type="outcome")
#> Warning: The following columns are not present but are helpful for harmonisation
#> eaf
pQTLtools::pqtlMR(exposure, outcome, prefix="IL6R-")
#> Harmonising IL.6 (O9YIMd) and Rheumatoid arthritis (9A6LZH)
#> Harmonising IL.6 (O9YIMd) and Coronary artery disease (GvLTWN)
#> Harmonising IL.6 (O9YIMd) and Atopic dermatitis (c4aS45)
#> Analysing 'O9YIMd' on '9A6LZH'
#> Analysing 'O9YIMd' on 'GvLTWN'
#> Analysing 'O9YIMd' on 'c4aS45'
pQTLtools::pqtlMR(exposure, outcome, prefix="IL6R_rev-",reverse=TRUE)
#> Harmonising IL.6 (O9YIMd) and Rheumatoid arthritis (9A6LZH)
#> Harmonising IL.6 (O9YIMd) and Coronary artery disease (GvLTWN)
#> Harmonising IL.6 (O9YIMd) and Atopic dermatitis (c4aS45)
#> Analysing '9A6LZH' on 'O9YIMd'
#> Analysing 'GvLTWN' on 'O9YIMd'
#> Analysing 'c4aS45' on 'O9YIMd'
unlink(c("IL6R*","pQTL-combined*"))
# Phenotype,SNP,effect_allele,other_allele,eaf,beta,se,pval
# ABO,rs505922,C,T,0.313,1.298,0.014,1.2e-1828
# LIFR,rs635634,T,C,0.180,-0.300,0.032,6.00E-21
# f <- file.path(find.package("pQTLtools",lib.loc=.libPaths()),"tests","ms.ins")
# exposure <- format_data(read.table(f, header=TRUE), samplesize_col="N")
# SNP Phenotype effect_allele other_allele eaf beta se pval N
# rs1800693 TNFB T C 0.6033 0.0282  0.0136 0.0389045   11787
# rs2364485 TNFB A C 0.1645 6514963 0.1759 1.62181e-20 11344
# https://raw.githubusercontent.com/MRCIEU/epigraphdb-pqtl/master/scripts/MR-pQTL-script.R
```
