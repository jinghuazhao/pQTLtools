# Limit of detection analysis

The function obtains lower limit of detection as in proteomic analysis.

## Usage

``` r
get.prop.below.LLOD(eset, flagged = "OUT")
```

## Arguments

-   eset:

    An ExpressionSet object.

-   flagged:

    A flag is an indicator for sample exclusion.

## Value

An updated ExpressionSet object.

## Author

James Peters

## Examples

``` r
suppressMessages(library(Biobase))
data(sample.ExpressionSet, package="Biobase")
exampleSet <- sample.ExpressionSet
Biobase::fData(exampleSet)
#> data frame with 0 columns and 500 rows
Biobase::fData(exampleSet)$lod.max <-
    apply(Biobase::exprs(exampleSet),1,quantile,runif(nrow(exampleSet)))
lod <- get.prop.below.LLOD(exampleSet)
x <- dplyr::arrange(fData(lod),desc(pc.belowLOD.new))
knitr::kable(head(lod))
#> 
#> 
#> |   | AFFX.MurIL2_at| AFFX.MurIL10_at| AFFX.MurIL4_at| AFFX.MurFAS_at| AFFX.BioB.5_at| AFFX.BioB.M_at|sex    |type    | score|
#> |:--|--------------:|---------------:|--------------:|--------------:|--------------:|--------------:|:------|:-------|-----:|
#> |A  |       192.7420|        97.13700|       45.81920|       22.54450|        96.7875|        89.0730|Female |Control |  0.75|
#> |B  |        85.7533|       126.19600|        8.83135|        3.60093|        30.4380|        25.8461|Male   |Case    |  0.40|
#> |C  |       176.7570|        77.92160|       33.06320|       14.68830|        46.1271|        57.2033|Male   |Control |  0.73|
#> |D  |       135.5750|        93.37130|       28.70720|       12.33970|        70.9319|        69.9766|Male   |Case    |  0.42|
#> |E  |        64.4939|        24.39860|        5.94492|       36.86630|        56.1744|        49.5822|Female |Case    |  0.93|
#> |F  |        76.3569|        85.50880|       28.29250|       11.25680|        42.6756|        26.1262|Male   |Control |  0.22|
#> |G  |       160.5050|        98.90860|       30.96940|       23.00340|        86.5156|        75.0083|Male   |Case    |  0.96|
#> |H  |        65.9631|        81.69320|       14.79230|       16.21340|        30.7927|        42.3352|Male   |Case    |  0.79|
#> |I  |        56.9039|        97.80150|       14.23990|       12.03750|        19.7183|        41.1207|Female |Case    |  0.37|
#> |J  |       135.6080|        90.48380|       34.48740|        4.54978|        46.3520|        91.5307|Male   |Control |  0.63|
#> |K  |        63.4432|        70.57330|       20.35210|        8.51782|        39.1326|        39.9136|Male   |Case    |  0.26|
#> |L  |        78.2126|        94.54180|       14.15540|       27.28520|        41.7698|        49.8397|Female |Control |  0.36|
#> |M  |        83.0943|        75.34550|       20.62510|       10.16160|        80.2197|        63.4794|Male   |Case    |  0.41|
#> |N  |        89.3372|        68.58270|       15.92310|       20.24880|        36.4903|        24.7007|Male   |Case    |  0.80|
#> |O  |        91.0615|        87.40500|       20.15790|       15.78490|        36.4021|        47.4641|Female |Case    |  0.10|
#> |P  |        95.9377|        84.45810|       27.81390|       14.32760|        35.3054|        47.3578|Female |Control |  0.41|
#> |Q  |       179.8450|        87.68060|       32.79110|       15.94880|        58.6239|        58.1331|Female |Case    |  0.16|
#> |R  |       152.4670|       108.03200|       33.52920|       14.67530|       114.0620|       104.1220|Male   |Control |  0.72|
#> |S  |       180.8340|       134.26300|       19.81720|       -7.91911|        93.4402|       115.8310|Male   |Case    |  0.17|
#> |T  |        85.4146|        91.40310|       20.41900|       12.88750|        22.5168|        58.1224|Female |Case    |  0.74|
#> |U  |       157.9890|        -8.68811|       26.87200|       11.91860|        48.6462|        73.4221|Male   |Control |  0.35|
#> |V  |       146.8000|        85.02120|       31.14880|       12.83240|        90.2215|        64.6066|Female |Control |  0.77|
#> |W  |        93.8829|        79.29980|       22.34200|       11.13900|        42.0053|        40.3068|Male   |Control |  0.27|
#> |X  |       103.8550|        71.65520|       19.01350|        7.55564|        57.5738|        41.8209|Male   |Control |  0.98|
#> |Y  |        64.4340|        64.23690|       12.16860|       19.98490|        44.8216|        46.1087|Female |Case    |  0.94|
#> |Z  |       175.6150|        78.70680|       17.37800|        8.96849|        61.7044|        49.4122|Female |Case    |  0.32|
plot(x[,2], main="Random quantile cut off", ylab="<lod%")
```
