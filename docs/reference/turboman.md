# Manhattan plots

Manhattan plots

## Usage

``` r
turboman(
  input_data_path,
  custom_peak_annotation_file_path,
  reference_file_path,
  pvalue_sign,
  plot_title,
  vertical_resolution = 1800
)
```

## Arguments

-   input_data_path:

    Path of the input association data.

-   custom_peak_annotation_file_path:

    Path of the custom annotation of variants.

-   reference_file_path:

    Path to the 'turboman_hg19_reference_data.rda' /
    'turboman_hg19_reference_data.rda' reference file.

-   pvalue_sign:

    Significance threshold p-value.

-   plot_title:

    Plot title which will be displayed on top of the plot.

-   vertical_resolution:

    A fixed number of points (pixels) to be plotted vertically.

## Details

**Input association data file / input_data_path**

Define the path of the input association data. The input data needs to
be a file that has:

1.  Spaces as field separators.

2.  One header line.

3.  Option I. (no extreme p-values present): 3 columns, being
    chromosome, position, pvalue in order, column names are not
    important. Option II. (extreme p-values present): 5 columns, being
    chromosome, position, pvalue, beta, se in order, column names are
    not important.

**Custom annotation file / custom_peak_annotation_file_path**

Define the path of the custom annotation of variants. The input data
needs to be a file that has:

1.  Spaces / tabs as field separators.

2.  One header line with exact column names (order not important).

3.  Columns: chromosome, position, label (e.g., gene name) /
    nearest_gene_name, cis/trans flag (optional). NB!: If no label is
    given, variants will be automatically annotated

**Reference file / reference_file_path**

Define the path to the 'turboman_hg19_reference_data.rda' /
'turboman_hg38_reference_data.rda' reference file that contains the LD
block breaks as in Berisa and Pickrell (2016) and gene coordinates used
to construct and annotate the Manhattan plot. Both are available from
the `turboman` directory of the installed package, e.g.,
`file.path(find.package('pQTLtools'),'turboman','turboman_hg38_reference_data.rda')`.

**Significance threshold p-value / pvalue_sign**

Define the significance threshold. This will be used to

1.  Highlight signal peaks that come above this significance threshold.

2.  Annotate the nearest gene to the top signal in the peak.

3.  Draw a horizontal reference line equal to this threshold.

**Title of the plot / plot_title**

Define title on top of the plot.

**Number of pixels on vertical axis / vertical_resolution**

Define a fixed number of points (pixels) on vertical axis.

## References

Berisa T, Pickrell JK (2016). “Approximately independent linkage
disequilibrium blocks in human populations.” *Bioinformatics*,
**32**(2), 283-285.
[doi:10.1093/bioinformatics/btv546](https://doi.org/10.1093/bioinformatics/btv546)
.

## Author

Arthur Gilly, Chris Finan, Bram Prins, see
<https://github.com/bpprins/turboman>.

## Examples

``` r
if (FALSE) { # \dontrun{
# Screen output
require(gap.datasets)
test <- mhtdata[c('chr','pos','p')]
write.table(test,file='test.txt',row.names=FALSE,quote=FALSE)
annotate <- subset(mhtdata[c('chr','start','gene','p')],p<5e-8 & gene!='')
names(annotate) <- c('chromosome','position','nearest_gene_name','p')
write.table(unique(annotate[,-4]),file='annotate.txt',row.names=FALSE,quote=FALSE)
input_data_path <- 'test.txt'
custom_peak_annotation_file_path <- 'annotate.txt'
reference_file_path <- file.path(find.package('pQTLtools',lib.loc=.libPaths()[1]),
  'turboman','turboman_hg19_reference_data.rda')
pvalue_sign <- 5e-8
plot_title <- 'gap.datasets example'
turboman(input_data_path, custom_peak_annotation_file_path,
         reference_file_path, pvalue_sign, plot_title)
# Figure shown on https://github.com/jinghuazhao/tests/tree/main/turboman
png('IL12B.png',width=3600, height=3600, pointsize = 12, res=450)
input_data_path <- 'IL.12B.txt.gz'
custom_peak_annotation_file_path <- 'IL.12B.annotate'
reference_file_path <-
  file.path(find.package('pQTLtools'),'turboman','turboman_hg19_reference_data.rda')
pvalue_sign <- 5e-8
plot_title <- 'IL12B'
turboman(input_data_path, custom_peak_annotation_file_path,
         reference_file_path, pvalue_sign, plot_title)
dev.off()
} # }
```
