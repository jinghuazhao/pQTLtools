# QQ plots

QQ plots

## Usage

``` r
turboqq(input_data_path, plot_title, plot_resolution = 1800)
```

## Arguments

-   input_data_path:

    Path of the input association data.

-   plot_title:

    Plot title to be displayed on top of the plot.

-   plot_resolution:

    a fixed number of points (pixels) to be plotted vertically and
    horizontally.

## Value

No direct return value. The script generates QQ plots as output.

## Details

The method uses the kth order statistic from a sample of n i.i.d. U(0,1)
statistics has a Beta(k,n+1-k) distribution as in Quesenberry and Hales
(1980) and Coded by Weale M, Price T.
<https://sites.google.com/site/mikeweale/software>

**Input association data file / input_data_path**

Define path of the input association data. The input data needs to be a
file that has:

1.  Spaces as field separators.

2.  One header line.

3.  Option I. (no extreme p-values present): 3 columns, being
    chromosome, position, pvalue in order, column names are not
    important. Option II. (extreme p-values present): 5 columns, being
    chromosome, position, pvalue, beta, se in order, column names are
    not important.

**Plot title / plot_title**

Define plot title which will be displayed on top of the plot.

**Resolution of plot / plot_resolution**

Define a fixed number of points (pixels) to be plotted vertically and
horizontally.

## References

Quesenberry CP, Hales C (1980). “Concentration bands for uniformity
plots.” *J Statist Comput Simul*, **11**(1), 41-53.
[doi:10.1080/00949658008810388](https://doi.org/10.1080/00949658008810388)
.

## Author

Bram Prins, <https://github.com/bpprins/turboqq>.

## Examples

``` r
if (FALSE) { # \dontrun{
png('test_qq.png', height = 1800, width = 1800, pointsize = 12, res = 450)
par(mar = c(4, 4, 3, 1))
require(gap.datasets)
test <- mhtdata[c('chr','pos','p')]
write.table(test,file='test.txt',row.names=FALSE,quote=FALSE)
input_data_path <- 'test.txt'
plot_title <- 'gap.datasets example'
turboqq(input_data_path, plot_title)
dev.off()
} # }
```
