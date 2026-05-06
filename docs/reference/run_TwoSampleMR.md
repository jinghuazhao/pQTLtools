# Basic TwoSampleMR analysis

Given harmonised data, this function conducts a two-sample MR analysis.

## Usage

``` r
run_TwoSampleMR(TwoSampleMRinput, mr_plot = "None", prefix = "")
```

## Arguments

-   TwoSampleMRinput:

    Harmonised data.

-   mr_plot:

    one of "None", "TwoSampleMR", "pQTLtools" for no, the original and
    the revised plots, respectively.

-   prefix:

    a prefix for output files.

## Value

No value is returned but several files.

## Details

As TwoSampleMR faces seemingly perplexing options, this function intends
to simplify various steps in a two-sample MR as in Dimou and Tsilidis
(2018) . It is particularly useful when a large numbher of MRs are
necessary, e.g., multiple proteins and their cis/trans regions need to
be examined, in which case prefix could direct the output to various
directories.

Check your authentication token if the example below fails to run.

## References

Dimou NL, Tsilidis KK (2018). “A Primer in Mendelian Randomization
Methodology with a Focus on Utilizing Published Summary Association
Data.” In Evangelou E (ed.), *Genetic Epidemiology: Methods and
Protocols*, chapter 13, 211–230. Springer New York, New York, NY. ISBN
978-1-4939-7868-7.
[doi:10.1007/978-1-4939-7868-7_13](https://doi.org/10.1007/978-1-4939-7868-7_13)
.

## Examples

``` r
suppressMessages(require(dplyr))
prot <- "MMP.10"
type <- "cis"
f <- paste0(prot,"-",type,".mrx")
d <- read.table(file.path(find.package("pQTLtools",lib.loc=.libPaths()),"tests",f),
                header=TRUE)
exposure <- TwoSampleMR::format_data(within(d,{P=10^logP}), phenotype_col="prot", snp_col="rsid",
                                     chr_col="Chromosome", pos_col="Posistion",
                                     effect_allele_col="Allele1", other_allele_col="Allele2",
                                     eaf_col="Freq1", beta_col="Effect", se_col="StdErr",
                                     pval_col="P", log_pval=FALSE,
                                     samplesize_col="N")
clump <- exposure[sample(1:nrow(exposure),nrow(exposure)/80),] # TwoSampleMR::clump_data(exposure)
# outcome <- TwoSampleMR::extract_outcome_data(snps=exposure$SNP,outcomes="ebi-a-GCST007432")
outcome <- pQTLtools::import_OpenGWAS("ebi-a-GCST007432","11:102090035-103364929","gwasvcf") %>%
           as.data.frame() %>%
           dplyr::mutate(outcome="FEV1",LP=10^-LP) %>%
           dplyr::select(ID,outcome,REF,ALT,AF,ES,SE,LP,SS,id) %>%
           setNames(c("SNP","outcome",paste0(c("other_allele","effect_allele","eaf","beta","se",
                                               "pval","samplesize","id"),".outcome")))
unlink("ebi-a-GCST007432.vcf.gz.tbi")
harmonise <- TwoSampleMR::harmonise_data(clump,outcome)
#> Harmonising MMP.10 (er44c0) and FEV1 (ebi-a-GCST007432)
prefix <- paste(prot,type,sep="-")
pQTLtools::run_TwoSampleMR(harmonise, mr_plot="pQTLtools", prefix=prefix)
#> Analysing 'er44c0' on 'ebi-a-GCST007432'
#> Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
#> ℹ Please use `linewidth` instead.
#> ℹ The deprecated feature was likely used in the pQTLtools package.
#>   Please report the issue at <https://github.com/jinghuazhao/pQTLtools/issues>.
#> Warning: The `size` argument of `element_line()` is deprecated as of ggplot2 3.4.0.
#> ℹ Please use the `linewidth` argument instead.
#> ℹ The deprecated feature was likely used in the pQTLtools package.
#>   Please report the issue at <https://github.com/jinghuazhao/pQTLtools/issues>.
#> `height` was translated to `width`.

#> `height` was translated to `width`.
#> Warning: Removed 1 row containing missing values or values outside the scale range
#> (`geom_point()`).


#> `height` was translated to `width`.
#> Warning: Removed 1 row containing missing values or values outside the scale range
#> (`geom_point()`).

caption <- "Table. MMP.10 variants and FEV1"
knitr::kable(read.delim(paste0(prefix,"-result.txt"),header=TRUE),
             caption=paste(caption, "(result)"))
#> 
#> 
#> Table: Table. MMP.10 variants and FEV1 (result)
#> 
#> |id.exposure |id.outcome       |outcome |exposure |method                    | nsnp|        b|      se|    pval|
#> |:-----------|:----------------|:-------|:--------|:-------------------------|----:|--------:|-------:|-------:|
#> |er44c0      |ebi-a-GCST007432 |FEV1    |MMP.10   |MR Egger                  |   10|  0.00526| 0.01677| 0.76191|
#> |er44c0      |ebi-a-GCST007432 |FEV1    |MMP.10   |Weighted median           |   10| -0.01824| 0.01056| 0.08403|
#> |er44c0      |ebi-a-GCST007432 |FEV1    |MMP.10   |Inverse variance weighted |   10| -0.02234| 0.00777| 0.00404|
#> |er44c0      |ebi-a-GCST007432 |FEV1    |MMP.10   |Simple mode               |   10| -0.05425| 0.03524| 0.15806|
#> |er44c0      |ebi-a-GCST007432 |FEV1    |MMP.10   |Weighted mode             |   10| -0.01338| 0.03957| 0.74306|
knitr::kable(read.delim(paste0(prefix,"-heterogeneity.txt"),header=TRUE),
             caption=paste(caption,"(heterogeneity)"))
#> 
#> 
#> Table: Table. MMP.10 variants and FEV1 (heterogeneity)
#> 
#> |id.exposure |id.outcome       |outcome |exposure |method                    |    Q| Q_df| Q_pval|
#> |:-----------|:----------------|:-------|:--------|:-------------------------|----:|----:|------:|
#> |er44c0      |ebi-a-GCST007432 |FEV1    |MMP.10   |MR Egger                  | 4.94|    8|  0.764|
#> |er44c0      |ebi-a-GCST007432 |FEV1    |MMP.10   |Inverse variance weighted | 8.39|    9|  0.496|
knitr::kable(read.delim(paste0(prefix,"-pleiotropy.txt"),header=TRUE),
             caption=paste(caption,"(pleiotropy)"))
#> 
#> 
#> Table: Table. MMP.10 variants and FEV1 (pleiotropy)
#> 
#> |id.exposure |id.outcome       |outcome |exposure | egger_intercept|      se| pval|
#> |:-----------|:----------------|:-------|:--------|---------------:|-------:|----:|
#> |er44c0      |ebi-a-GCST007432 |FEV1    |MMP.10   |        -0.00368| 0.00198|  0.1|
knitr::kable(read.delim(paste0(prefix,"-single.txt"),header=TRUE),
             caption=paste(caption,"(single)"))
#> 
#> 
#> Table: Table. MMP.10 variants and FEV1 (single)
#> 
#> |exposure |outcome |id.exposure |id.outcome       | samplesize|SNP                             |        b|      se|       p|
#> |:--------|:-------|:-----------|:----------------|----------:|:-------------------------------|--------:|-------:|-------:|
#> |MMP.10   |FEV1    |er44c0      |ebi-a-GCST007432 |     321047|rs11601744                      | -0.01770| 0.02056| 0.38918|
#> |MMP.10   |FEV1    |er44c0      |ebi-a-GCST007432 |     321047|rs140957819                     | -0.03994| 0.03391| 0.23889|
#> |MMP.10   |FEV1    |er44c0      |ebi-a-GCST007432 |     321047|rs2292733                       |  0.00228| 0.01901| 0.90448|
#> |MMP.10   |FEV1    |er44c0      |ebi-a-GCST007432 |     321047|rs2514414                       | -0.05914| 0.03226| 0.06675|
#> |MMP.10   |FEV1    |er44c0      |ebi-a-GCST007432 |     321047|rs484316                        | -0.05898| 0.03155| 0.06154|
#> |MMP.10   |FEV1    |er44c0      |ebi-a-GCST007432 |     321047|rs577850                        |  0.00870| 0.03285| 0.79124|
#> |MMP.10   |FEV1    |er44c0      |ebi-a-GCST007432 |     321047|rs608419                        | -0.06043| 0.03309| 0.06784|
#> |MMP.10   |FEV1    |er44c0      |ebi-a-GCST007432 |     321047|rs6590984                       | -0.01855| 0.01408| 0.18781|
#> |MMP.10   |FEV1    |er44c0      |ebi-a-GCST007432 |     321047|rs75524134                      |  0.00284| 0.03354| 0.93246|
#> |MMP.10   |FEV1    |er44c0      |ebi-a-GCST007432 |     321047|rs9735884                       | -0.05507| 0.03478| 0.11335|
#> |MMP.10   |FEV1    |er44c0      |ebi-a-GCST007432 |     321047|All - Inverse variance weighted | -0.02234| 0.00777| 0.00404|
#> |MMP.10   |FEV1    |er44c0      |ebi-a-GCST007432 |     321047|All - MR Egger                  |  0.00526| 0.01677| 0.76191|
knitr::kable(read.delim(paste0(prefix,"-loo.txt"),header=TRUE),
             caption=paste(caption,"(loo)"))
#> 
#> 
#> Table: Table. MMP.10 variants and FEV1 (loo)
#> 
#> |exposure |outcome |id.exposure |id.outcome       | samplesize|SNP         |       b|      se|       p|
#> |:--------|:-------|:-----------|:----------------|----------:|:-----------|-------:|-------:|-------:|
#> |MMP.10   |FEV1    |er44c0      |ebi-a-GCST007432 |     321047|rs11601744  | -0.0231| 0.00856| 0.00695|
#> |MMP.10   |FEV1    |er44c0      |ebi-a-GCST007432 |     321047|rs140957819 | -0.0214| 0.00803| 0.00783|
#> |MMP.10   |FEV1    |er44c0      |ebi-a-GCST007432 |     321047|rs2292733   | -0.0273| 0.00851| 0.00136|
#> |MMP.10   |FEV1    |er44c0      |ebi-a-GCST007432 |     321047|rs2514414   | -0.0201| 0.00801| 0.01217|
#> |MMP.10   |FEV1    |er44c0      |ebi-a-GCST007432 |     321047|rs484316    | -0.0200| 0.00802| 0.01273|
#> |MMP.10   |FEV1    |er44c0      |ebi-a-GCST007432 |     321047|rs577850    | -0.0242| 0.00800| 0.00250|
#> |MMP.10   |FEV1    |er44c0      |ebi-a-GCST007432 |     321047|rs608419    | -0.0201| 0.00799| 0.01185|
#> |MMP.10   |FEV1    |er44c0      |ebi-a-GCST007432 |     321047|rs6590984   | -0.0240| 0.00948| 0.01136|
#> |MMP.10   |FEV1    |er44c0      |ebi-a-GCST007432 |     321047|rs75524134  | -0.0238| 0.00799| 0.00292|
#> |MMP.10   |FEV1    |er44c0      |ebi-a-GCST007432 |     321047|rs9735884   | -0.0206| 0.00797| 0.00969|
#> |MMP.10   |FEV1    |er44c0      |ebi-a-GCST007432 |     321047|All         | -0.0223| 0.00777| 0.00404|
for (x in c("result","heterogeneity","pleiotropy","single","loo"))
    unlink(paste0(prefix,"-",x,".txt"))
```
