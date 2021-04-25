# pQTL MR based on MR-Base

File `Instruments.csv` was converted from its `Instruments.xlsx` and can be read with
`Ins <- read.csv("Instruments.csv"),Phenotype %in% c("ABO","LIFR"))` and a simplified
version is obtained with `awk -FS="," 'NR==1||$1~/ABO/||$1~/LIFR/' ~/R/work/Instruments.csv > Ins.csv`.

```r
Ins <- read.csv("Ins.csv")
ids <- c("ieu-a-7","ieu-a-798","ukb-a-61")
pqtlMR(Ins,ids)
```

The original script by Jie Zheng
```r
library("readxl")
Ins <- data <- read_excel("Instruments.xlsx",1)
ids <- scan("outcome.id.txt",what="")
```

* [Instruments.xlsx](https://github.com/MRCIEU/epigraphdb-pqtl/blob/master/data/Instruments.xlsx?raw=true)
* [outcome.id.txt](https://raw.githubusercontent.com/MRCIEU/epigraphdb-pqtl/master/data/outcome.id.txt)
* [MR-pQTL-script.R](https://raw.githubusercontent.com/MRCIEU/epigraphdb-pqtl/master/scripts/MR-pQTL-script.R)

## A test

```r
abc <- function(a=1,b=2,c=3) list(sum=a+b+c,prod=a*b*c)
test <- function(x,y=x$sum,z=x$prod,u=y+z) u
test(abc())
```
