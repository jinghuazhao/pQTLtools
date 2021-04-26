# pQTL MR (instrumental variables) and outcomes (MR-Base)

## Experiment

File `Instruments.csv` was converted from its `Instruments.xlsx` and can be read with
`Ins <- read.csv("Instruments.csv"),Phenotype %in% c("ABO","LIFR"))` and or directly
from `awk -FS="," 'NR==1||$1~/ABO/||$1~/LIFR/' Instruments.csv > Ins.csv`.

```r
source("pqtlMR.R")
Ins <- read.csv("Ins.csv")
ids <- c("ieu-a-7","ebi-a-GCST007432")
pqtlMR(Ins,ids)
```

The original script by Jie Zheng

```r
library("readxl")
Ins <- data <- read_excel("Instruments.xlsx",1)
ids <- scan("outcome.id.txt",what="")
```
handles three files

1. [Instruments.xlsx](https://github.com/MRCIEU/epigraphdb-pqtl/blob/master/data/Instruments.xlsx?raw=true)
2. [outcome.id.txt](https://raw.githubusercontent.com/MRCIEU/epigraphdb-pqtl/master/data/outcome.id.txt)
3. [MR-pQTL-script.R](https://raw.githubusercontent.com/MRCIEU/epigraphdb-pqtl/master/scripts/MR-pQTL-script.R)

## A toy example

This offers possibility to simplify two-sample MR,

```r
abc <- function(a=1,b=2,c=3) list(sum=a+b+c,prod=a*b*c)
test <- function(x,y=x$sum,z=x$prod,u=y+z) u
test(abc())
```
but unfortunately the scheme does not always work.

---

## Legacy implementation

These can be seen from `utils.R` and pqtlMR.Rd and `run_TwoSampleMR.Rd`.
