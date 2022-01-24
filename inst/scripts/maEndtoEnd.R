settup <- function()
{
## ----biocSetup, include = FALSE, results="hide", warning=FALSE, echo=FALSE----
## Changing the YAML to the following changes the output to 'latex'
## output:
## BiocWorkflowTools::f1000_article
## More at:
## https://stackoverflow.com/questions/35144130/in-knitr-how-can-i-test-for-if-the-output-will-be-pdf-or-word

## add figures to current directory if Rmd is converted to LaTeX
## This simplifies processing on overleaf

if (!is.null(knitr::opts_knit$get("rmarkdown.pandoc.to"))) {
    
    to_html <- knitr::opts_knit$get("rmarkdown.pandoc.to") != 'latex'
        if (!to_html) {
          knitr::opts_chunk$set(fig.path = "", out.width = '.8\\linewidth')
        } 

}
## here, pdf specific adapations can be made
## however, the option "rmarkdown.pandoc.to" is only set by rmarkdown::render,
## so it returns NULL if run interactively!


## ----testAlternativeToChunkAbove, eval=FALSE, include=FALSE, echo=FALSE-------
#  # this does not work as it will only test the output formats
#  # specified in the YAML header
#  any(grepl("f1000", rmarkdown::all_output_formats(knitr::current_input())))

## ----knitrOptions, include=FALSE----------------------------------------------
library(BiocStyle)
library(knitr)
options(digits = 4, width = 80, timeout = 600)
opts_chunk$set(echo = TRUE, 
               tidy = FALSE, 
               include = TRUE,
               dev = c('png'),
               fig.width = 6, fig.height = 3.5,
               comment = '  ', 
               dpi = 300, cache = TRUE,
               fig.pos = "H")

}
## ----installBioc, eval=FALSE--------------------------------------------------
#  if (!require("BiocManager"))
#      install.packages("BiocManager")
#  BiocManager::install("maEndToEnd", version = "devel")

## ---- echo=TRUE, results="hide", warning=FALSE, eval=FALSE--------------------
#  
#  #In order to download the package from GitHub, we need the "install_github"
#  #function from the "remotes" package. We download the latest developer
#  #version of "remotes" from GitHub with the devtool::install_github
#  #function; note that this is necessary as the current "remotes" version on
#  #CRAN doesn't allow us to correctly download the "maEndToEnd" package:
#  
#  install.packages("devtools")
#  library(devtools)
#  
#  devtools::install_github("r-lib/remotes")
#  library(remotes)
#  packageVersion("remotes") # has to be 1.1.1.9000 or later
#  
#  remotes::install_github("b-klaus/maEndToEnd", ref="master")
#  

## ----maEndToEndImport---------------------------------------------------------
suppressPackageStartupMessages({library("maEndToEnd")})

## ----pkgList, results="hide"--------------------------------------------------
#General Bioconductor packages
    library(Biobase)
    library(oligoClasses)
     
#Annotation and data import packages
    library(ArrayExpress)
    library(pd.hugene.1.0.st.v1)
    library(hugene10sttranscriptcluster.db)
     
#Quality control and pre-processing packages
    library(oligo)
    library(arrayQualityMetrics)
     
#Analysis and statistics packages
    library(limma)
    library(topGO)
    library(ReactomePA)
    library(clusterProfiler)
     
#Plotting and color options packages
    library(gplots)
    library(ggplot2)
    library(geneplotter)
    library(RColorBrewer)
    library(pheatmap)
    library(enrichplot)
     
#Formatting/documentation packages
   #library(rmarkdown)
   #library(BiocStyle)
    library(dplyr)
    library(tidyr)

#Helpers:
    library(stringr)
    library(matrixStats)
    library(genefilter)
    library(openxlsx)
   #library(devtools)


## ----generateFolderForRawData-------------------------------------------------
raw_data_dir <- tempdir()

if (!dir.exists(raw_data_dir)) {
    dir.create(raw_data_dir)
}


## ----getDataEBI, eval=TRUE, results='hide', message=FALSE---------------------
anno_AE <- getAE("E-MTAB-2967", path = raw_data_dir, type = "raw")

## ----sumexp, fig.cap = "Structure of Bioconductor’s ExpressionSet class.", echo=FALSE, fig.show="asis", fig.width = 7, fig.height = 4.5----
par(mar=c(0,0,0,0))
plot(1,1,xlim=c(0,100),ylim=c(0, 150),bty="n",
     type="n",xlab="",ylab="",xaxt="n",yaxt="n")

polygon(c(50,85,85,50),c(70,70,130,130),col="lightcoral",border=NA)
polygon(c(55, 56, 56, 55),c(70,70,130,130),col=rgb(1,0,0,.5),border=NA)
polygon(c(50,85,85,50),c(118,118,120,120),col=rgb(1,0,0,.5),border=NA)
text(67.5,100,"assay(s)", cex = 1)
text(67.5,90,"e.g. 'exprs'", cex = 1)

polygon(c(45,49,49,45),c(70,70,130,130),col="honeydew2",border=NA)
text(47,100, "microarray probes", srt=90, cex=1)

polygon(c(50,85,85,50),c(131,131,140,140),col="darkseagreen3",border=NA)
text(67.5,135.5,"sample IDs", cex = 1)


polygon(c(20,40,40,20),c(70,70,130,130),col=rgb(0,0,1,.5),border=NA)
polygon(c(20,40,40,20),c(118,118,120,120),col=rgb(0,0,1,.5),border=NA)
text(30,100,"featureData", cex = 1)

polygon(c(15,19,19,15),c(70,70,130,130),col="honeydew2",border=NA)
text(17,100, "microarray probes", srt=90, cex=1)

polygon(c(20,40,40,20),c(131,131,140,140),col="bisque3",border=NA)
text(30,135.5,"features", cex = 1)

polygon(c(50,65, 65, 50),c(55, 55, 0, 0),col=rgb(.5,0,.5,.5),border=NA)
polygon(c(55, 56, 56, 55),c(55, 55, 0, 0),col=rgb(.5,0,.5,.5),border=NA)
text(57.5, 30,"phenoData", cex = 1, srt=270)

polygon(c(50,65, 65, 50),c(56,56,65,65),col="darkseagreen3",border=NA)
text(57.5,60.5,"sample IDs", cex = 1)

## ----getSDRF------------------------------------------------------------------
sdrf_location <- file.path(raw_data_dir, "E-MTAB-2967.sdrf.txt")
SDRF <- read.delim(sdrf_location)

rownames(SDRF) <- SDRF$Array.Data.File
SDRF <- AnnotatedDataFrame(SDRF)

## ----importCelfiles, results="hide", eval=TRUE, dependson="getSDRF", warning = FALSE----
raw_data <- oligo::read.celfiles(filenames = file.path(raw_data_dir, 
                                                SDRF$Array.Data.File),
                                    verbose = FALSE, phenoData = SDRF)
stopifnot(validObject(raw_data))

## ----inspectPhenoData, eval=TRUE----------------------------------------------
head(Biobase::pData(raw_data))

## ----reassignpData, eval=TRUE-------------------------------------------------
Biobase::pData(raw_data) <- Biobase::pData(raw_data)[, c("Source.Name",
                                     "Characteristics.individual.",
                                     "Factor.Value.disease.",
                                     "Factor.Value.phenotype.")]

## ----inspectAssayData, eval=TRUE----------------------------------------------
Biobase::exprs(raw_data)[1:5, 1:5]

## ----qualityControlRawDataPCA, fig.cap="PCA plot of the log–transformed raw expression data."----
exp_raw <- log2(Biobase::exprs(raw_data))
PCA_raw <- prcomp(t(exp_raw), scale. = FALSE)

percentVar <- round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dataGG <- data.frame(PC1 = PCA_raw$x[,1], PC2 = PCA_raw$x[,2],
                    Disease = pData(raw_data)$Factor.Value.disease.,
                    Phenotype = pData(raw_data)$Factor.Value.phenotype.,
                    Individual = pData(raw_data)$Characteristics.individual.)

ggplot(dataGG, aes(PC1, PC2)) +
      geom_point(aes(shape = Disease, colour = Phenotype)) +
  ggtitle("PCA plot of the log-transformed raw expression data") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5))+
  coord_fixed(ratio = sd_ratio) +
  scale_shape_manual(values = c(4,15)) + 
  scale_color_manual(values = c("darkorange2", "dodgerblue4"))

## ----qualityControlRawDataBox, fig.cap="Intensity boxplots of the log2–transformed raw data."----
oligo::boxplot(raw_data, target = "core", 
               main = "Boxplot of log2-intensitites for the raw data")

## ----arrayQualityMetricsRaw, eval = TRUE, warning=FALSE, message=FALSE--------
arrayQualityMetrics(expressionset = raw_data,
    outdir = tempdir(),
    force = TRUE, do.logtransform = TRUE,
    intgroup = c("Factor.Value.disease.", "Factor.Value.phenotype."))

## ----annotationDataBaseContent, eval = TRUE-----------------------------------
head(ls("package:hugene10sttranscriptcluster.db"))

## ----DifferenceBetweenExonAndGeneTypeArrays, fig.width=7, fig.height=6, echo=FALSE, fig.cap="Visualization of the difference between \"Exon\" type array (left) and \"Gene\" type array (right)."----
library(grid)

par(mar=c(0,0,0,0))
plot(1,1,xlim=c(0,100),ylim=c(0, 20),bty="n",
     type="n",xlab="",ylab="",xaxt="n",yaxt="n")
dat <- data.frame(x = rep(seq(0, 2.9, 0.1), 30), 
                  y = rep(seq(0, 2.9, 0.1), each = 30))


col <- c("blue", "darkslategray1", "goldenrod1", "yellow", "royalblue1", "coral")
dat$colours<-sample(rep(col, 150))
dat$colours2 <- sample(c(rep("yellow", 2), rep("goldenrod1", 5), rep("coral",7), 
                  rep("blue", 5), rep ("darkslategray1",3), rep("royalblue1", 8), rep(gray(0:9/10), 87)))

# opening the graphic device and
# setting up a viewport with borders:
vp1 <- viewport(x = 0.1, y = 0.1, w = 0.12, h = 0.12, 
                just = c("left", "bottom"), name = "vp1")


# plotting rectangles using x/y positions
grid.rect(x=dat$x,y=dat$y, height=0.1, width=0.1, hjust=0,vjust=0,vp=vp1,
          gp=gpar(col=1, fill=as.character(dat$colours)))

##########################

vp2 <- viewport(x = 0.5, y = 0.1, w = 0.12, h = 0.12,
                just = c("left", "bottom"), name = "vp1")


# plotting rectangles using x/y positions
grid.rect(x=dat$x,y=dat$y, height=0.1, width=0.1, hjust=0,vjust=0,vp=vp2,
          gp=gpar(col=1, fill=as.character(dat$colours2)))

## ----RMAcalibrationForRLE, eval=TRUE------------------------------------------
palmieri_eset <- oligo::rma(raw_data, target = "core", normalize = FALSE)

## ----boxplotDataForRLETidy, fig.cap="Boxplot for the RLE values", warning=FALSE----
row_medians_assayData <- 
  Biobase::rowMedians(as.matrix(Biobase::exprs(palmieri_eset)))

RLE_data <- sweep(Biobase::exprs(palmieri_eset), 1, row_medians_assayData)

RLE_data <- as.data.frame(RLE_data)
RLE_data_gathered <- 
  tidyr::gather(RLE_data, patient_array, log2_expression_deviation)

ggplot2::ggplot(RLE_data_gathered, aes(patient_array,
                                       log2_expression_deviation)) + 
  geom_boxplot(outlier.shape = NA) + 
  ylim(c(-2, 2)) + 
  theme(axis.text.x = element_text(colour = "aquamarine4", 
                                  angle = 60, size = 6.5, hjust = 1 ,
                                  face = "bold"))

## ----RMAcalibrationWITHnormalization, eval=TRUE-------------------------------
palmieri_eset_norm <- oligo::rma(raw_data, target = "core")

## ----PCAMetricsCalibrated, fig.cap = "PCA plot of the calibrated, summarized data.", eval = TRUE----
exp_palmieri <- Biobase::exprs(palmieri_eset_norm)
PCA <- prcomp(t(exp_palmieri), scale = FALSE)

percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2],
                    Disease = 
                     Biobase::pData(palmieri_eset_norm)$Factor.Value.disease.,
                    Phenotype = 
                     Biobase::pData(palmieri_eset_norm)$Factor.Value.phenotype.)


ggplot(dataGG, aes(PC1, PC2)) +
      geom_point(aes(shape = Disease, colour = Phenotype)) +
  ggtitle("PCA plot of the calibrated, summarized data") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_fixed(ratio = sd_ratio) +
  scale_shape_manual(values = c(4,15)) + 
  scale_color_manual(values = c("darkorange2", "dodgerblue4"))

## ----rownamesForHeatmap, fig.height = 8.5, fig.width = 7, eval = TRUE, echo=TRUE----
phenotype_names <- ifelse(str_detect(pData
                                    (palmieri_eset_norm)$Factor.Value.phenotype.,
                             "non"), "non_infl.", "infl.")

disease_names <- ifelse(str_detect(pData
                                    (palmieri_eset_norm)$Factor.Value.disease.,
                             "Crohn"), "CD", "UC")

annotation_for_heatmap <- 
  data.frame(Phenotype = phenotype_names,  Disease = disease_names)

row.names(annotation_for_heatmap) <- row.names(pData(palmieri_eset_norm))


## ----HeatmapWithAnnotation, fig.height = 8.5, fig.width = 7, fig.cap="Heatmap of the sample-to-sample distances"----
dists <- as.matrix(dist(t(exp_palmieri), method = "manhattan"))

rownames(dists) <- row.names(pData(palmieri_eset_norm))
hmcol <- rev(colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"))(255))
colnames(dists) <- NULL
diag(dists) <- NA

ann_colors <- list(
  Phenotype = c(non_infl. = "chartreuse4", infl. = "burlywood3"),
  Disease = c(CD = "blue4", UC = "cadetblue2")
                   )
pheatmap(dists, col = (hmcol), 
         annotation_row = annotation_for_heatmap,
         annotation_colors = ann_colors,
         legend = TRUE, 
         treeheight_row = 0,
         legend_breaks = c(min(dists, na.rm = TRUE), 
                         max(dists, na.rm = TRUE)), 
         legend_labels = (c("small distance", "large distance")),
         main = "Clustering heatmap for the calibrated samples")

## ----IntensityBasedManualFiltering, fig.cap="Histogram of the median intensities per gene"----
palmieri_medians <- rowMedians(Biobase::exprs(palmieri_eset_norm))

hist_res <- hist(palmieri_medians, 100, col = "cornsilk1", freq = FALSE, 
            main = "Histogram of the median intensities", 
            border = "antiquewhite4",
            xlab = "Median intensities")

## ----setManualThreshold, fig.cap="Histogram of the median intensities per gene with manual intensity filtering threshold (red line)."----
man_threshold <- 4

hist_res <- hist(palmieri_medians, 100, col = "cornsilk", freq = FALSE, 
            main = "Histogram of the median intensities",
            border = "antiquewhite4",
            xlab = "Median intensities")

abline(v = man_threshold, col = "coral4", lwd = 2)

## ----expGroups, dependson="PCAMetricsCalibrated"------------------------------
no_of_samples <- 
  table(paste0(pData(palmieri_eset_norm)$Factor.Value.disease., "_", 
                  pData(palmieri_eset_norm)$Factor.Value.phenotype.))
no_of_samples 

## ----filteringOfLowIntensity_transcripts--------------------------------------
samples_cutoff <- min(no_of_samples)

idx_man_threshold <- apply(Biobase::exprs(palmieri_eset_norm), 1,
                           function(x){
                          sum(x > man_threshold) >= samples_cutoff})
                          table(idx_man_threshold)

palmieri_manfiltered <- subset(palmieri_eset_norm, idx_man_threshold)

## ----annotateData, eval=TRUE, dependson="intensityBasedFiltering", message = FALSE----
anno_palmieri <- AnnotationDbi::select(hugene10sttranscriptcluster.db,
                                  keys = (featureNames(palmieri_manfiltered)),
                                  columns = c("SYMBOL", "GENENAME"),
                                  keytype = "PROBEID")

anno_palmieri <- subset(anno_palmieri, !is.na(SYMBOL))

## ----multipleMappings, dependson="annotateData"-------------------------------
anno_grouped <- group_by(anno_palmieri, PROBEID)
anno_summarized <- 
  dplyr::summarize(anno_grouped, no_of_matches = n_distinct(SYMBOL))

head(anno_summarized)

anno_filtered <- filter(anno_summarized, no_of_matches > 1)

head(anno_filtered)

probe_stats <- anno_filtered 

nrow(probe_stats)

## ----excludeMultipleMappingsFromAssayData, dependson="multipleMappings"-------
ids_to_exlude <- (featureNames(palmieri_manfiltered) %in% probe_stats$PROBEID)

table(ids_to_exlude)

palmieri_final <- subset(palmieri_manfiltered, !ids_to_exlude)

validObject(palmieri_final)

## ----recallAnnoPalmieri-------------------------------------------------------
head(anno_palmieri)

## ----excludeMultipleMappingsFromFeatureData-----------------------------------
fData(palmieri_final)$PROBEID <- rownames(fData(palmieri_final))

## ----excludeMultipleMappingsFromFeatureData2----------------------------------
fData(palmieri_final) <- left_join(fData(palmieri_final), anno_palmieri)


# restore rownames after left_join
rownames(fData(palmieri_final)) <- fData(palmieri_final)$PROBEID 
    
validObject(palmieri_final)

## ----introductionOfAbbreviations----------------------------------------------
individual <- 
  as.character(Biobase::pData(palmieri_final)$Characteristics.individual.)

tissue <- str_replace_all(Biobase::pData(palmieri_final)$Factor.Value.phenotype.,
                  " ", "_")

tissue <- ifelse(tissue == "non-inflamed_colonic_mucosa",
                 "nI", "I")

disease <- 
  str_replace_all(Biobase::pData(palmieri_final)$Factor.Value.disease.,
                  " ", "_")
disease <- 
  ifelse(str_detect(Biobase::pData(palmieri_final)$Factor.Value.disease., 
                    "Crohn"), "CD", "UC")

## ----createDesign, eval=TRUE, dependson="excludeMultipleMappings"-------------
i_CD <- individual[disease == "CD"]
design_palmieri_CD <- model.matrix(~ 0 + tissue[disease == "CD"] + i_CD)
colnames(design_palmieri_CD)[1:2] <- c("I", "nI")
rownames(design_palmieri_CD) <- i_CD 

i_UC <- individual[disease == "UC"]
design_palmieri_UC <- model.matrix(~ 0 + tissue[disease == "UC"] + i_UC )
colnames(design_palmieri_UC)[1:2] <- c("I", "nI")
rownames(design_palmieri_UC) <- i_UC 

## ----inspectDesignMatrix, eval = TRUE, dependson="createDesign"---------------
head(design_palmieri_CD[, 1:6])

head(design_palmieri_UC[, 1:6])

## ----visualizeExpressionChanges, fig.cap="Visualization of expression changes"----
tissue_CD <- tissue[disease == "CD"]
crat_expr <- Biobase::exprs(palmieri_final)["8164535", disease == "CD"]
crat_data <- as.data.frame(crat_expr)
colnames(crat_data)[1] <- "org_value"
crat_data <- mutate(crat_data, individual = i_CD, tissue_CD)

crat_data$tissue_CD <- factor(crat_data$tissue_CD, levels = c("nI", "I"))

ggplot(data = crat_data, aes(x = tissue_CD, y = org_value, 
                             group = individual, color = individual)) +
      geom_line() +
      ggtitle("Expression changes for the CRAT gene")

## ----obtainFitCRAT------------------------------------------------------------
crat_coef <- lmFit(palmieri_final[,disease == "CD"],
                design = design_palmieri_CD)$coefficients["8164535",]

crat_coef

## ----2obtainFitCRAT-----------------------------------------------------------
crat_fitted <- design_palmieri_CD %*% crat_coef
rownames(crat_fitted) <- names(crat_expr)
colnames(crat_fitted) <- "fitted_value"

crat_fitted

## ----3obtainFitCRAT, fig.cap="Expression changes for the CRAT gene"-----------
crat_data$fitted_value <- crat_fitted

ggplot(data = crat_data, aes(x = tissue_CD, y = fitted_value, 
                             group = individual, color = individual)) +
      geom_line() +
      ggtitle("Fitted expression changes for the CRAT gene")

## ----tTest--------------------------------------------------------------------
crat_noninflamed <- na.exclude(crat_data$org_value[tissue == "nI"])
crat_inflamed <- na.exclude(crat_data$org_value[tissue == "I"])
res_t <- t.test(crat_noninflamed ,crat_inflamed , paired = TRUE)
res_t

## ----createContrastMatrixAndFitModel, eval=TRUE, dependson="createDesign"-----
contrast_matrix_CD <- makeContrasts(I-nI, levels = design_palmieri_CD)

palmieri_fit_CD <- eBayes(contrasts.fit(lmFit(palmieri_final[,disease == "CD"],
                                design = design_palmieri_CD),
                                contrast_matrix_CD))

contrast_matrix_UC <- makeContrasts(I-nI, levels = design_palmieri_UC)

palmieri_fit_UC <- eBayes(contrasts.fit(lmFit(palmieri_final[,disease == "UC"],
                                design = design_palmieri_UC),
                                contrast_matrix_UC))

## ----extractResultsCD, eval = TRUE, dependson="createContrastMatrixAndFitModel", message=FALSE, fig.cap="Histogram of the p–values for Crohn’s disease."----
table_CD <- topTable(palmieri_fit_CD, number = Inf)
head(table_CD)

hist(table_CD$P.Value, col = brewer.pal(3, name = "Set2")[1],
     main = "inflamed vs non-inflamed - Crohn's disease", xlab = "p-values")


## ----extractResultsUC, eval = TRUE, dependson="createContrastMatrixAndFitModel", message=FALSE, fig.cap="Histogram of the p–values for ulcerative colitis."----

table_UC <- topTable(palmieri_fit_UC, number = Inf)
head(table_UC)

hist(table_UC$P.Value, col = brewer.pal(3, name = "Set2")[2],
     main = "inflamed vs non-inflamed - Ulcerative colitis", xlab = "p-values")


## ----pCutForCD----------------------------------------------------------------
nrow(subset(table_UC, P.Value < 0.001))

## ----FDRforUC-----------------------------------------------------------------
tail(subset(table_UC, P.Value < 0.001))

## ----compareDEgenes, dependson=c("extractResultsUC", "extractResultsCD")------
fpath <- system.file("extdata", "palmieri_DE_res.xlsx", package = "maEndToEnd")
palmieri_DE_res <- sapply(1:4, function(i) read.xlsx(cols = 1, fpath, 
                                                     sheet = i, startRow = 4))

names(palmieri_DE_res) <- c("CD_UP", "CD_DOWN", "UC_UP", "UC_DOWN")
palmieri_DE_res <- lapply(palmieri_DE_res, as.character)
paper_DE_genes_CD <- Reduce("c", palmieri_DE_res[1:2])
paper_DE_genes_UC <- Reduce("c", palmieri_DE_res[3:4])

overlap_CD <- length(intersect(subset(table_CD, P.Value < 0.001)$SYMBOL,  
                               paper_DE_genes_CD)) / length(paper_DE_genes_CD)


overlap_UC <- length(intersect(subset(table_UC, P.Value < 0.001)$SYMBOL,
                               paper_DE_genes_UC)) / length(paper_DE_genes_UC)
overlap_CD
overlap_UC 

total_genenumber_CD <- length(subset(table_CD, P.Value < 0.001)$SYMBOL)
total_genenumber_UC <- length(subset(table_UC, P.Value < 0.001)$SYMBOL)

total_genenumber_CD
total_genenumber_UC

## ----VolcanoPlot, fig.cap="Volcano plot of the DE-genes", fig.height=8, fig.width=7----

volcano_names <- ifelse(abs(palmieri_fit_CD$coefficients)>=1, 
                        palmieri_fit_CD$genes$SYMBOL, NA)
             
             
volcanoplot(palmieri_fit_CD, coef = 1L, style = "p-value", highlight = 100, 
            names = volcano_names,
            xlab = "Log2 Fold Change", ylab = NULL, pch=16, cex=0.35)

## ----FDRcontrolledDEgenes, dependson=c("extractResultsUC", "extractResultsCD"), eval=TRUE----
DE_genes_CD <- subset(table_CD, adj.P.Val < 0.1)$PROBEID

## ----GOAnalysisCreateBackgrounds, eval=TRUE, warning=FALSE, message=FALSE-----
back_genes_idx <- genefilter::genefinder(palmieri_final, 
                                        as.character(DE_genes_CD), 
                                        method = "manhattan", scale = "none")

## ----GOAnalysisCreateBackgrounds2, eval=TRUE, warning=FALSE, message=FALSE----
back_genes_idx <- sapply(back_genes_idx, function(x)x$indices)

## ----GOAnalysisCreateBackgrounds3, eval=TRUE, warning=FALSE, message=FALSE----
back_genes <- featureNames(palmieri_final)[back_genes_idx]
back_genes <- setdiff(back_genes, DE_genes_CD)

    
intersect(back_genes, DE_genes_CD)
length(back_genes)

## ----multidensityPlot, fig.cap="Selecting a background set of genes for the gene ontology analysis."----
multidensity(list(
        all = table_CD[,"AveExpr"] ,
        fore = table_CD[DE_genes_CD , "AveExpr"],
        back = table_CD[rownames(table_CD) %in% back_genes, "AveExpr"]),
        col = c("#e46981", "#ae7ee2", "#a7ad4a"),
     xlab = "mean expression",
   main = "DE genes for CD-background-matching")

## ----createFactorOfInterestingGenes, dependson="GOAnalysisCreateBackgrounds", eval=TRUE----
gene_IDs <- rownames(table_CD)
in_universe <- gene_IDs %in% c(DE_genes_CD, back_genes)
in_selection <- gene_IDs %in% DE_genes_CD 

all_genes <- in_selection[in_universe]
all_genes <- factor(as.integer(in_selection[in_universe]))
names(all_genes) <- gene_IDs[in_universe] 

## ----createTopGODataSet, dependson="createFactorOfInterestingGenes", eval=TRUE, message = FALSE, warning=FALSE----
top_GO_data <- new("topGOdata", ontology = "BP", allGenes = all_genes,
 nodeSize = 10, annot = annFUN.db, affyLib = "hugene10sttranscriptcluster.db")

## ----runtopGOTests, results='hide', eval=TRUE, dependson = "createTopGODataSet", message = FALSE----
result_top_GO_elim <- 
  runTest(top_GO_data, algorithm = "elim", statistic = "Fisher")
result_top_GO_classic <- 
  runTest(top_GO_data, algorithm = "classic", statistic = "Fisher")

## ----processtopGOResults, eval=TRUE, dependson="runtopGOTests"----------------
res_top_GO <- GenTable(top_GO_data, Fisher.elim = result_top_GO_elim,
        Fisher.classic = result_top_GO_classic,
        orderBy = "Fisher.elim" , topNodes = 100)

genes_top_GO <- printGenes(top_GO_data, whichTerms = res_top_GO$GO.ID,
    chip = "hugene10sttranscriptcluster.db", geneCutOff = 1000)

res_top_GO$sig_genes <- sapply(genes_top_GO, function(x){
                str_c(paste0(x[x$'raw p-value' == 2, "Symbol.id"],";"), 
                      collapse = "")
    })

head(res_top_GO[,1:8], 20)

## ----graphOfResults, fig.height = 6, eval=TRUE, results='hide', dpi=600, fig.cap="Significantly enriched GO nodes in the GO hierarchy"----
# showSigOfNodes(top_GO_data, score(result_top_GO_elim), firstSigNodes = 3,
#                useInfo = 'def')

## ----mapIDsToEntrez, dependson="createFactorOfInterestingGenes", message = FALSE----
entrez_ids <- mapIds(hugene10sttranscriptcluster.db, 
      keys = rownames(table_CD), 
      keytype = "PROBEID",
      column = "ENTREZID")

## ----runReactomeEnrichment, fig.cap="Enriched Reactome pathways and their p–values as a bar chart.", eval = TRUE, warning=FALSE----
reactome_enrich <- enrichPathway(gene = entrez_ids[DE_genes_CD], 
                                universe = entrez_ids[c(DE_genes_CD, 
                                                        back_genes)],
                                organism = "human",
                                pvalueCutoff = 0.05,
                                qvalueCutoff = 0.9, 
                                readable = TRUE)

reactome_enrich@result$Description <- paste0(str_sub(
                                    reactome_enrich@result$Description, 1, 20),
                                    "...")

head(as.data.frame(reactome_enrich))[1:6]

## ----reactomeBar, dependson="runReactomeEnrichment", fig.cap="Enriched Reactome pathways and their p–values as a bar chart."----
barplot(reactome_enrich)

## ----emapplot, dependson="runReactomeEnrichment", fig.width=6, fig.height = 7, fig.cap="Enriched Reactome pathways enrichment results as a graph."----
reactome_enrich <- pairwise_termsim(reactome_enrich)
emapplot(reactome_enrich, showCategory = 10)

## -----------------------------------------------------------------------------
gc()

length(getLoadedDLLs())

sessionInfo()
