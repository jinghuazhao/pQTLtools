#' Manhattan plots
#'
#' @param input_data_path Path of the input association data
#' @param output_data_rootname Root name of the plot output file
#' @param custom_peak_annotation_file_path Path of the custom annotation of variants
#' @param reference_file_path Path to the "turboman_hg19_reference_data.rda" reference file
#' @param pvalue_sign Significance threshold p-value
#' @param plot_title Plot title which will be displayed on top of the plot
#' @export
#' @return NULL
#'
#' @details
#' **Input association data path / input_data_path**
#'
#' Define the path of the input association data...
#' The input data needs to be a file that has:
#' 1. Spaces as field separators
#' 2. One header line
#' 3. Two option:
#'     I. (no extreme p-values present): 3 columns, being chromosome, position, pvalue in order, column names are not important.
#'     II. (extreme p-values present): 5 columns, being chromosome, position, pvalue, beta, se in order, column names are not important.
#'
#' **Output data root name / output_data_rootname**
#'
#' Define the root name of the plot output file...
#' Ex.: "my_man_plot" will result in an output file named "my_man_plot.png"
#'
#' **OPTIONAL: Custom annotation file / custom_peak_annotation_file_path**
#'
#' Define the path of the custom annotation of variants...
#' The input data needs to be a file that has:
#' 1. Spaces / tabs as field separators
#' 2. One header line with exact column names (order not important)
#' 3. 3 columns: chromosome, position, label (e.g., gene name)
#'    Column names chromosome, position, nearest_gene_name
#' NB!: If no label is given, variants will be automatically annotated
#'
#' **Reference file path / reference_file_path**
#'
#' Define the path to the "turboman_hg19_reference_data.rda" reference file that
#' contains the LD block breaks as in \insertCite{berisa16;textual}{pQTLtools}
#' and gene coordinates used to construct and annotate the Manhattan plot
#'
#' **Significance threshold p-value / pvalue_sign**
#'
#' Define the significance threshold...
#' This will be used to...
#' 1. Highlight signal peaks that come above this significance threshold
#' 2. Annotate the nearest gene to the top signal in the peak
#' 3. Draw a horizontal reference line equal to this threshold
#'
#' @examples
#' \dontrun{
#' require(gap.datasets)
#' test <- mhtdata[c("chr","pos","p")]
#' write.table(test,file="test.txt",row.names=FALSE,quote=FALSE)
#' annotate <- unique(subset(mhtdata[c("chr","start","gene","p")],p<5e-8 & gene!=""))
#' names(annotate) <- c("chromosome","position","nearest_gene_name","p")
#' write.table(annotate[,-4],file="annotate.txt",row.names=FALSE,quote=FALSE)
#' input_data_path <- "test.txt"
#' output_data_rootname <- "test_man"
#' custom_peak_annotation_file_path <- "annotate.txt"
#' reference_file_path <- "turboman_hg19_reference_data.rda"
#' pvalue_sign <- 5e-8
#' plot_title <- "gap.datasets example"
#' turboman(input_data_path, output_data_rootname, custom_peak_annotation_file_path, reference_file_path, pvalue_sign, plot_title)
#'
#' input_data_path <- "FCN3.txt.gz"
#' output_data_rootname <- "FCN3"
#' custom_peak_annotation_file_path <- "FCN3.txt"
#' reference_file_path <- "turboman_hg19_reference_data.rda"
#' pvalue_sign <- 5e-8
#' plot_title <- "FCN3"
#' turboman(input_data_path, output_data_rootname, custom_peak_annotation_file_path, reference_file_path, pvalue_sign, plot_title)
#' }
#'
#' @author Arthur Gilly, Chris Finan, Bram Prins
#'
#' **Changelog**
#' * 2017/12/04. Added Varying window size in which gene names are plotted based on the length of gene names so that gene names will not be truncated in the plot
#' * 2017/12/04. Added an option allowing the user to supply a table with top SNPs to highlight, which can, for example, be selected by conditional analyses / LD pruning
#' * 2017/12/04. Fixed a bug (found by Tao) that leads to erroneous binning when plotting extremely sparse data
#' * 2017/12/05. Added optional log-pval calculations from beta/SE for extreme p-values that run below `.Machine$double.xmin`
#' * 2017/12/05. Added an automated option for plotting already log-transformed p-values (user can supply either)
#' * 2017/12/05. Increased contrast between chromosomes on the plot (Eric H)
#' * 2017/12/05. Set '10' in the y-axis label (-log10 pvalue) to subscript (Eric H)
#'
#' **To do list**
#' 1. Display gene names bold / non-bold when genic / intergenic
#' 2. Create LD blocks for hg38
#' 3. Add many more sanity checks (are you giving pvalues, the file does not exist, etc.)
#' 4. Create a function version of the script (Eric H)
#' 5. ...
#'
#' @references
#' \insertAllCited{}

turboman <- function(input_data_path, output_data_rootname, custom_peak_annotation_file_path, reference_file_path, pvalue_sign, plot_title) {
###================================================================================================================================================###
###	1. Defining input settings                                                                                                                 ###
###================================================================================================================================================###

## Verbose printing status bars
fat_status_bar<-strrep("=",108)
skinny_status_bar<-strrep("-",108)
print_status_bar <- function(message) {cat("\n", fat_status_bar, "\n", message, "\n", fat_status_bar, "\n\n", sep = "")}
print_status_bar("1. Defining input settings")
## Log start time
start.time <- Sys.time()
print(paste0(" Starting at ",start.time))
## Read in arguments from the command line
for (arg in commandArgs(trailingOnly=TRUE))
{
  ta = strsplit(arg,"=",fixed=TRUE)
  if(!is.na(ta[[1]][2]))
  {
    assign(ta[[1]][1],ta[[1]][2])
  } else {
    stop("Not all arguments are given")
  }
}
## Check if the custom peak annotation data file exists 
custom_peak_annotation_file_path_exists<-exists("custom_peak_annotation_file_path")
## Assign variable classes
input_data_path <- as.character(input_data_path)
output_data_rootname <- as.character(output_data_rootname)
if (custom_peak_annotation_file_path_exists) {
custom_peak_annotation_file_path <- as.character(custom_peak_annotation_file_path)
}
reference_file_path <- as.character(reference_file_path)
pvalue_sign <- as.numeric(pvalue_sign)
plot_title <- as.character(plot_title)
## Load the reference data
load(reference_file_path)

###================================================================================================================================================###
### 2. Reading in the plotting data, log10 transform the pvalues, initial basic data sanity checks                                                 ###
###================================================================================================================================================###

print_status_bar("2. Reading in the plotting data, log10 transform the pvalues, initial basic data sanity checks")
print(paste0("  Data file path : ",input_data_path))

#---------------------------------------------------------------------
# Reading in association plot data with scan
#---------------------------------------------------------------------

initial_data_dims<-dim(as.data.frame(read.table(input_data_path, header=TRUE, stringsAsFactors=FALSE, nrows=10)))[2]
if (initial_data_dims==3) {
initial_data <- data.frame(scan(input_data_path,
                        what = list(chromosome = 0, position = 0,pvalue= 0),
                        skip=1,
                        sep=" ",
                        quiet=TRUE))
initial_data_contains_beta_se<-FALSE
} else if (initial_data_dims==5) {
initial_data <- data.frame(scan(input_data_path,
                        what = list(chromosome = 0, position = 0,pvalue= 0, beta=0, se=0),
                        skip=1,
                        sep=" ",
                        quiet=TRUE))
initial_data_contains_beta_se<-TRUE
} else {
  stop("Input data does not have expected dimensions")
  }

#---------------------------------------------------------------------
# Reading in peak annotation plot data with read.table and check if 
# it is already annotated with labels / gene names (# Check if the 
# variants are annotated with gene names
#---------------------------------------------------------------------

if (custom_peak_annotation_file_path_exists) {
   gene_plot_data <- data.frame(read.table(custom_peak_annotation_file_path,header=TRUE,stringsAsFactors=FALSE))
   nearest_gene_names_annotated<-("nearest_gene_name" %in% colnames(gene_plot_data))
} else {
   nearest_gene_names_annotated<-FALSE
}

#---------------------------------------------------------------------
# Preparing the association data
#---------------------------------------------------------------------

## Check if p-values are already logged
if (length(which(initial_data$pvalue>1))>0) {
    initial_data$log_pvalue<-initial_data$pvalue
    ## Get only the complete data
    initial_data<-initial_data[complete.cases(initial_data),]
    # Remove the original pvalues
    initial_data$pvalue<-NULL
} else {
    ## Calculate the -log10 p-value for the input data 
    initial_data$log_pvalue<--log10(initial_data$pvalue)
    ## If beta/SE are provided, and pvalues are missing (because they are extreme), log10 P recalculate from beta/SE
    missing_pvalues_index<-which((is.na(initial_data$log_pvalue) | initial_data$log_pvalue==0))
    if (initial_data_contains_beta_se & (length(missing_pvalues_index)>0)) {
        # Calculate expected p-values for missing data
        missing_pvalues<-(-log(2, base=10)-pnorm(-abs(initial_data$beta[missing_pvalues_index]/initial_data$se[missing_pvalues_index]), log=T)/log(10))
        # Only replace if indeed they were below the smallest non-zero normalized floating-point number
        initial_data[missing_pvalues_index,c("log_pvalue")]<-ifelse(missing_pvalues > -log10(.Machine$double.xmin),missing_pvalues,NA)
    }
    ## Get only the complete data
    initial_data<-initial_data[complete.cases(initial_data),]
    # Remove the original pvalues
    initial_data$pvalue<-NULL
}

# Calculate -log10 of significance threshold
log_pvalue_sign<--log10(pvalue_sign)

###================================================================================================================================================###
### 3. Preparing data for plotting, calculating variables related to plotting                                                                      ###
###================================================================================================================================================###

print_status_bar("3. Preparing data for plotting, calculating variables related to plotting")
## Set vertical resolution
# We will maximally allow a fixed number of points to be plotted vertically, 
# here choosing 800 as a 'pixel' unit on a normal standard R plot.
vertical_resolution<-1800
## Find the largest p-value, which we will use to make the y-axis 'resolution'
observed_log_pvalue_maximum<-max(initial_data$log_pvalue,na.rm=TRUE)
## Now we will scale the 800-point resolution for the p-values 
log_pvalue_break_size<-observed_log_pvalue_maximum/vertical_resolution
## Create a vector from 0 to the vertical resolution, which we will use to bin pvalues
scaling_vector<-seq(0,vertical_resolution,by=log_pvalue_break_size)
## Obtain all unique chromosome numbers for which we have data
chromosomes<-unlist(unique(initial_data$chromosome),recursive = FALSE, use.names = FALSE)

###================================================================================================================================================###
### 4. Data reduction procedure, creating plotting data                                                                                            ###
###================================================================================================================================================###

print_status_bar("4. Data reduction procedure, creating plotting data")
## Create an empty dataframe for the plotting data
plot_data<-NULL
plot_data<-as.data.frame(plot_data)
## If the custom annotation data is not given, create an empty data frame for annotation downstream which will be filled
if (!custom_peak_annotation_file_path_exists) { 
  gene_plot_data<-NULL
  gene_plot_data<-as.data.frame(gene_plot_data)
}

## Initialise a counter that keeps track of how many top signals we have
top_snp_counter=0
## Define the LD block data (loaded from rda reference data file)
ld_block_breaks<-ld_block_breaks_pickrell_hg19_eur
## Define the gene annotation data (loaded from rda reference data file)
gene_coordinates<-refgene_gene_coordinates_h19
## Start the loop that will go over all unique chromosomes to reduce the data
for (chromosome_number in chromosomes) {
## Verbose progress tracker
johnny_bravo <- ifelse (chromosome_number%%2==0, " ha !", " hoo !")
print (paste0("chromosome ",chromosome_number,johnny_bravo))

#---------------------------------------------------------------------
# Reduce the position data for plotting using the LD blocks (X-axis)
#---------------------------------------------------------------------

## Extract the LD block breaks for the chromosome
chromosome_ld_block_breaks<-ld_block_breaks[which(ld_block_breaks[,1]==chromosome_number),2]
## Count the number of LD block breaks for the chromosome
number_of_ld_block_bins<-as.integer(length(chromosome_ld_block_breaks)-1)
## Select the association data for this chromosome to be reduced
initial_data_chromosome<-initial_data[which(initial_data$chromosome==chromosome_number),]
## Assign bin numbers to the positions based on the LD block breaks
initial_data_chromosome$bin<-findInterval(initial_data_chromosome$position,chromosome_ld_block_breaks)
## Calculate the midpoints of the breaks, which we will use as X-coordinates on the Manhattan plot.
plot_x_coordinates<-(head(chromosome_ld_block_breaks, -1) + diff(chromosome_ld_block_breaks) / 2)
## Create a temporary dataframe in which we will first assemble the plotting data per chromosome
plot_data_per_chromosome_df<-NULL
plot_data_per_chromosome_df<-as.data.frame(plot_data_per_chromosome_df)
## Determine the number of bins in a chromosome, for which we will loop over to reduce the p_value data
unique_x_bins<-number_of_ld_block_bins
for (bin_number in 1:unique_x_bins){
 
#-----------------------------------------------------------------------------------------
# Reduce the p-value data for plotting to imagined resolution of max 800 points vertically
#-----------------------------------------------------------------------------------------
    
    ## Create a temporary dataframe in which we will first assemble the plotting data per bin in a chromosome
    plot_data_per_bin_in_chromosome_df<-NULL
    plot_data_per_bin_in_chromosome_df<-as.data.frame(plot_data_per_bin_in_chromosome_df)
    ## Now reduce (bin) the p-values (Y-axis values) to 800 bins, and multiply each bin (starting from 1 to max resolution) by the calculated pvalue_break_size
    plot_data_per_bin_in_chromosome_pvalues<-(unique(.bincode((initial_data_chromosome[['log_pvalue']][which(initial_data_chromosome$bin==bin_number)]), scaling_vector, right = TRUE, include.lowest = FALSE)*log_pvalue_break_size))
    ## If there are no p-values for a bin, enter one line with chromosome and position, missing pvalue, and missing highlight value
    if(length(plot_data_per_bin_in_chromosome_pvalues)==0){
    plot_data_per_bin_in_chromosome_df[1,1]<-chromosome_number
    plot_data_per_bin_in_chromosome_df[1,2]<-plot_x_coordinates[bin_number]
    plot_data_per_bin_in_chromosome_df[1,3]<-NA
    plot_data_per_bin_in_chromosome_df[1,4]<-NA
    ## Name the columns in the plotting data dataframe that was made for this bin
    colnames(plot_data_per_bin_in_chromosome_df)<-c("chromosome","position","log_pvalue","highlight_vector")
    } else {
    ## First find the top SNP in the bin based on the maximum p-value in this bin
    largest_pvalue_in_bin<-max((initial_data_chromosome[['log_pvalue']][which(initial_data_chromosome$bin==bin_number)]),na.rm=TRUE)[1]
    largest_pvalue_index_in_bin<-which(initial_data_chromosome[['log_pvalue']][which(initial_data_chromosome$bin==bin_number)]==largest_pvalue_in_bin)
    ## If the largest p-value in the bin is significant, continue with finding the matching position for the SNP with highest p-value, and it's nearest gene
    ## but only if a custom annotation file is not given
    if ((largest_pvalue_in_bin>log_pvalue_sign) & (!custom_peak_annotation_file_path_exists)) {
      ## Increase top SNP counter
      top_snp_counter=top_snp_counter+1
      ## Define the position of the top SNP in the bin
      largest_pvalue_index_in_bin_position<-initial_data_chromosome[which(initial_data_chromosome$bin==bin_number),][largest_pvalue_index_in_bin,][['position']][1]
      ## Extract the gene annotation data from the gene table for this particular chromosome
      gene_coordinates_chromosome<-gene_coordinates[which(gene_coordinates$chromosome==chromosome_number),]
      ## Find the smallest distances to the position of our top SNP
      #OLD CODE : smallest_distance_to_gene_for_top_snp_in_bin<-min(abs(gene_coordinates_chromosome$gene_transcription_midposition-largest_pvalue_index_in_bin_position),na.rm=TRUE)
      smallest_distance_to_gene_start_for_top_snp_in_bin<-min(abs(gene_coordinates_chromosome$gene_transcription_start-largest_pvalue_index_in_bin_position),na.rm=TRUE)
      smallest_distance_to_gene_stop_for_top_snp_in_bin<-min(abs(gene_coordinates_chromosome$gene_transcription_stop-largest_pvalue_index_in_bin_position),na.rm=TRUE)
      ## Find which gene corresponds to the smallest distances to the position of our top SNP
      # OLD CODE : genename_for_top_snp_in_bin<-as.character(gene_coordinates_chromosome[which(abs(gene_coordinates_chromosome$gene_transcription_midposition-largest_pvalue_index_in_bin_position)==smallest_distance_to_gene_for_top_snp_in_bin),c("gene_name")])[1]
      if (smallest_distance_to_gene_start_for_top_snp_in_bin<smallest_distance_to_gene_stop_for_top_snp_in_bin) {
            genename_for_top_snp_in_bin<-as.character(gene_coordinates_chromosome[which(abs(gene_coordinates_chromosome$gene_transcription_start-largest_pvalue_index_in_bin_position)==smallest_distance_to_gene_start_for_top_snp_in_bin),c("gene_name")])[1]
      } else {
            genename_for_top_snp_in_bin<-as.character(gene_coordinates_chromosome[which(abs(gene_coordinates_chromosome$gene_transcription_stop-largest_pvalue_index_in_bin_position)==smallest_distance_to_gene_stop_for_top_snp_in_bin),c("gene_name")])[1]
      }
      ## Enter the chromosome of the top SNP in the gene annotation dataframe which we will use in the plot
      gene_plot_data[top_snp_counter,1]<-chromosome_number
      ##  Enter the mid-bin coordinate for the top SNP in the gene annotation dataframe which we will use in the plot
      gene_plot_data[top_snp_counter,2]<-plot_x_coordinates[bin_number]
      ## Enter the pvalue for the top SNP in the gene annotation dataframe which we will use in the plot
      gene_plot_data[top_snp_counter,3]<-largest_pvalue_in_bin
      ## Enter the nearest gene for the top SNP in the gene annotation dataframe which we will use in the plot
      gene_plot_data[top_snp_counter,4]<-genename_for_top_snp_in_bin
      ## Name the columns of the gene annotation dataframe
      colnames(gene_plot_data)<-c("chromosome","position","log_pvalue","nearest_gene_name")
      }
    ## Else if there are p-values for a bin, enter the chromosome, the midposition for this bin, the pvalue, and vector values telling whether this bin
    ## should be highlighted in the plot 
    plot_data_per_bin_in_chromosome_df[1:length(plot_data_per_bin_in_chromosome_pvalues),1]<-rep(chromosome_number,length(plot_data_per_bin_in_chromosome_pvalues))
    plot_data_per_bin_in_chromosome_df[1:length(plot_data_per_bin_in_chromosome_pvalues),2]<-rep(plot_x_coordinates[bin_number],length(plot_data_per_bin_in_chromosome_pvalues))
    plot_data_per_bin_in_chromosome_df[1:length(plot_data_per_bin_in_chromosome_pvalues),3]<-plot_data_per_bin_in_chromosome_pvalues
    plot_data_per_bin_in_chromosome_df[1:length(plot_data_per_bin_in_chromosome_pvalues),4]<-ifelse(largest_pvalue_in_bin>log_pvalue_sign,
                                                                                                    rep(1,length(plot_data_per_bin_in_chromosome_pvalues)), 
                                                                                                    rep(0,length(plot_data_per_bin_in_chromosome_pvalues)))
    ## Name the columns in the plotting data dataframe that was made for this bin
    colnames(plot_data_per_bin_in_chromosome_df)<-c("chromosome","position","log_pvalue","highlight_vector")
    }
    ## Add the per-bin plotting data dataframe to the per-chromosome plotting data dataframe
    plot_data_per_chromosome_df<-rbind(plot_data_per_chromosome_df,plot_data_per_bin_in_chromosome_df)
    ## Remove the per-bin plotting data dataframe
    rm(plot_data_per_bin_in_chromosome_df)
}
## Add the per-chromosome plotting data dataframe to the per-chromosome plotting data dataframe
plot_data<-rbind(plot_data,plot_data_per_chromosome_df)
## Remove the per-chromosome plotting data dataframe
rm(plot_data_per_chromosome_df)
rm(initial_data_chromosome)
}

#-----------------------------------------------------------------------------------------
# Processing custom annotation for variants 
#-----------------------------------------------------------------------------------------

## If the custom peak annotation file exists (chromosomes and positions of variants provided),
## but there are no labels present, perform an annotation for nearest genes and look up their
## pvalues in the association data

if  (custom_peak_annotation_file_path_exists & !nearest_gene_names_annotated) {
    ## Create an additional column filled with NA
    gene_plot_data$nearest_gene_name<-NA
    ## Count how many variants should be annotated
    number_of_peak_annotations<-dim(gene_plot_data)[1]
    ## Loop over the variants, find the nearest genes and p-values of these variants in the association data
    for (peak_number in 1:number_of_peak_annotations){
        peak_snp_chromosome<-gene_plot_data$chromosome[peak_number]
        peak_snp_position<-gene_plot_data$position[peak_number]
        ## Extract the gene annotation data from the gene table for this particular chromosome
        gene_coordinates_chromosome<-gene_coordinates[which(gene_coordinates$chromosome==peak_snp_chromosome),]
        ## Find the smallest distances to the position of our top SNP
        smallest_distance_to_gene_for_top_snp_in_bin<-min(abs(gene_coordinates_chromosome$gene_transcription_midposition-peak_snp_position),na.rm=TRUE)
        ## Find which gene corresponds to the smallest distances to the position of our top SNP
        genename_for_top_snp_in_bin<-as.character(gene_coordinates_chromosome[which(abs(gene_coordinates_chromosome$gene_transcription_midposition-peak_snp_position)==smallest_distance_to_gene_for_top_snp_in_bin),c("gene_name")])[1]
        gene_plot_data$nearest_gene_name[peak_number]<-genename_for_top_snp_in_bin
        ## Also find the p-vals !!!
        gene_plot_data$log_pvalue[peak_number]<-initial_data[which(initial_data$chromosome==peak_snp_chromosome & initial_data$position==peak_snp_position),c("log_pvalue")][1]
    }
} else if (custom_peak_annotation_file_path_exists & nearest_gene_names_annotated) {
    ## Count how many variants should be annotated
    number_of_peak_annotations<-dim(gene_plot_data)[1]
    ## Loop over the variants, find the nearest genes and p-values of these variants in the association data
    for (peak_number in 1:number_of_peak_annotations){
        peak_snp_chromosome<-gene_plot_data$chromosome[peak_number]
        peak_snp_position<-gene_plot_data$position[peak_number]
        ## Also find the p-vals !!!
        gene_plot_data$log_pvalue[peak_number]<-initial_data[which(initial_data$chromosome==peak_snp_chromosome & initial_data$position==peak_snp_position),c("log_pvalue")][1]
    }
} else {
    cat("\n Custom annotation data fully provided\n\n")
}

###================================================================================================================================================###
### 5. Start plotting                                                                                                                             ###
###================================================================================================================================================###

print_status_bar("5. Start plotting")
## Define image properties
png(paste0(output_data_rootname,".png"),width=3600, height=3600, pointsize = 12, res=450)
# pdf(paste0(output_data_rootname,".pdf"), height=5, width=8)
## Sorting the plotting data from the GWAS datafile and the gene annotation file
plot_data<-plot_data[order(plot_data$chromosome,plot_data$position),]
if (dim(gene_plot_data)[1] > 0 ) {
    gene_plot_data<-gene_plot_data[order(gene_plot_data$chromosome,gene_plot_data$position),]
}
## Putting the data in vectors, easier to work with
chromosomes <- plot_data$chromosome
positions <- plot_data$position
log_pvalues <- plot_data$log_pvalue
highlight_vector <- plot_data$highlight_vector
unique_chromosomes<-unique(chromosomes)
## Defining the Y-axis maxima to be used for the Y-axis limit to plot the association data and the gene annotation
## Truncate the maximum p-value to an integer
log_pvalue_truncated <- trunc(max(log_pvalues,na.rm=TRUE))
# Use the truncated maximum p-value to define nice Y-axis limits
y_axis_plot_data_limit<-ifelse(log_pvalue_truncated < 10,10, 
                        ifelse(log_pvalue_truncated < 15, 15,
                        ifelse(log_pvalue_truncated < 15, 15,
                        (((trunc(log_pvalue_truncated/10)+1)*10)))))
## Define the limits for the gene annotations / lines and start of gene name display and real size of plot window
y_axis_stop_gene_annotation_vertical_lines<-y_axis_plot_data_limit
y_axis_stop_gene_annotation_diagonal_lines<-y_axis_plot_data_limit*1.1
y_axis_true_limit<-y_axis_plot_data_limit*1.3
### Setting the plotting data
## How many chromosomes / plot only the chromosomes for which there's data
max_nchr<-length(unique(chromosomes))
## Set indices
x <- 1:max_nchr
x2<- 1:max_nchr
## Define the actual plotting X coordinates as will be used in the plot, based on basepair positions of the variants
for (i in 1:max_nchr)
{
     chromosome_number=which(chromosomes==i)
     x[i] <- trunc((max(na.omit(positions[chromosome_number])))/100) +100000
     x2[i] <- trunc((min(na.omit(positions[chromosome_number])))/100) -100000
}
x[1]=x[1]-x2[1]
x2[1]=0-x2[1]
for (i in 2:(max_nchr+1))
{
	x[i] <- x[i-1]-x2[i]+x[i]
	x2[i] <- x[i-1]-x2[i]
}
## Calculate the final x-coordinates of the association data to plot
x_coordinates = trunc(positions/100) + x2[chromosomes]
## Define the x-axis limit
x_axis_limit <-max(x_coordinates,na.rm=TRUE)-min(x_coordinates,na.rm=TRUE)
## Set the final y-coordinates of the association data to plot
y_coordinates = log_pvalues
## Calculate the final x-coordinates of the top SNPs with annotated genes to plot
if (dim(gene_plot_data)[1] > 0 ) {
    gene_x_coordinates = trunc(gene_plot_data$position/100) + x2[gene_plot_data$chromosome]
    ## Set the final y-coordinates of the top SNPs with annotated genes to plot
    gene_y_coordinates = gene_plot_data$log_pvalue
    ## Set the nearest gene names of the top SNPs to plot
    nearest_gene_names_hits = gene_plot_data$nearest_gene_name
    nearest_gene_names_cistrans = gene_plot_data$cistrans
}
## Set the default colors of the all association datapoints, grey and light grey, alternating between odd and even chromosome numbers
col1="gray72"
col2="gray50"
chromosome_colour <- ifelse (chromosomes%%2==0, col1, col2)
## Plot the association data
plot(x_coordinates,y_coordinates,pch=20,col=chromosome_colour,axes=F,ylab="",xlab="",bty="n",ylim=c(0,y_axis_true_limit),cex=0.8,main=plot_title)
## Plot the gene annotation data
if (dim(gene_plot_data)[1] > 0 ) {
    ## Create X-axis breaks at which the gene names will be plotted, chosing random number of 150 as I simply
    ## Assume maximally 150 peaks will be annotated and at 150 genes I hope no gene names will be displayed overlapping
    x_axis_break_factor<-x_axis_limit/150
    ## Draw top SNP gene annotation lines
    y_axis_stop_gene_annotation_vertical_lines<-y_axis_plot_data_limit
    y_axis_stop_gene_annotation_diagonal_lines<-y_axis_plot_data_limit*1.1
    y_axis_true_limit<-y_axis_plot_data_limit*1.3
    for (i in 1:length(gene_x_coordinates))
    {
    ## Define the coordinates for the vertical annotation lines
    vertical_annotation_line_x_coordinate_start<-gene_x_coordinates[i]
    vertical_annotation_line_x_coordinate_stop<-gene_x_coordinates[i]
    vertical_annotation_line_y_coordinate_start<-gene_y_coordinates[i]
    vertical_annotation_line_y_coordinate_stop<-y_axis_stop_gene_annotation_vertical_lines
    ## Define the diagonal for the vertical annotation lines
    diagonal_annotation_line_x_coordinate_start<-gene_x_coordinates[i]
    diagonal_annotation_line_x_coordinate_stop<-((x_axis_limit/length(gene_x_coordinates))/2)+((i-1)*(x_axis_limit/length(gene_x_coordinates)))
    diagonal_annotation_line_y_coordinate_start<-y_axis_stop_gene_annotation_vertical_lines
    diagonal_annotation_line_y_coordinate_stop<-y_axis_stop_gene_annotation_diagonal_lines
    ## Draw the vertical annotation lines
    lines(c(vertical_annotation_line_x_coordinate_start,vertical_annotation_line_x_coordinate_stop),c(vertical_annotation_line_y_coordinate_start,vertical_annotation_line_y_coordinate_stop),col="grey", lty=2,lwd=1)
    ## Draw the diagonal annotation lines
    lines(c(diagonal_annotation_line_x_coordinate_start,diagonal_annotation_line_x_coordinate_stop),c(diagonal_annotation_line_y_coordinate_start,diagonal_annotation_line_y_coordinate_stop),col="grey", lty=2,lwd=1)
    ## Plot the gene names for each top SNP
    # calculating font sizes 
    number_of_annotations_to_plot<-dim(gene_plot_data)[1]
    maximum_characters_annotation<-max(nchar(gene_plot_data$nearest_gene_name),na.rm=TRUE)
    if (( number_of_annotations_to_plot <= 70) & (maximum_characters_annotation <=9)) {
       gene_label_cex_size<-1
    } else if (( number_of_annotations_to_plot > 70) & (maximum_characters_annotation <=9)) {
       gene_label_cex_size<-1.30-(0.006*number_of_annotations_to_plot)
    } else if (( number_of_annotations_to_plot <= 70) & (maximum_characters_annotation >9)) {
       gene_label_cex_size<-1.15-(0.03264*maximum_characters_annotation)
    } else (( number_of_annotations_to_plot > 70) & (maximum_characters_annotation >9))
       gene_label_cex_size_n_genes_annotation<-1.30-(0.006*number_of_annotations_to_plot)
       gene_label_cex_size_max_char_annotation<-1.15-(0.03264*maximum_characters_annotation)
       gene_label_cex_size<-ifelse(gene_label_cex_size_n_genes_annotation<gene_label_cex_size_max_char_annotation,
                                   gene_label_cex_size_n_genes_annotation,gene_label_cex_size_max_char_annotation)
    # Plot the labels
    text(diagonal_annotation_line_x_coordinate_stop,diagonal_annotation_line_y_coordinate_stop,labels=nearest_gene_names_hits[i],cex=gene_label_cex_size,srt=90,adj = c(0,0.5),font=3,ps=12,col=ifelse(nearest_gene_names_cistrans[i]=="cis","red","blue"))
    }
}
## Draw the significanc threshold
lines(c(0,max(x_coordinates,na.rm=TRUE)),c(log_pvalue_sign,log_pvalue_sign),col="dodgerblue4", lty=2,lwd=1)
## Draw horizontal axis(can't find a way to nicely draw a x-axis without not running over points :/
lines(c(0,max(x_coordinates,na.rm=TRUE)),c((y_axis_plot_data_limit/200)*-1,(y_axis_plot_data_limit/200)*-1),col="black", lty=1,lwd=1)
## Highlight the significant peaks
points(x_coordinates[which(highlight_vector==1)], y_coordinates[which(highlight_vector==1)], col="dodgerblue4", pch=20, cex=0.8)
## Make (non-existing) x-axis text and labels
x_axis_chromosome_labels_text<-c(1:22,"X","Y XY M","","")
x_axis_chromosome_labels<-x_axis_chromosome_labels_text[unique_chromosomes]
for (i in 1:max_nchr)
{
    label_positions = (x[i] + x2[i]) / 2
    mtext(x_axis_chromosome_labels[i],1,at=label_positions,cex=1,line=0,las=2)
}
## Set the horizontal axis text - "Chromosome" and each chromosome number for which the data is plotted
mtext("Chromosome",1,at=x[max_nchr]/2,cex=1,line=2)
## Draw the y-axis with value ticks
axis(2,las=2,pos=0,yaxp=c(0,y_axis_plot_data_limit,10))
## Draw the y-axis label
mtext(expression(paste(-"log"[10], " p-value")),2,line=1)
## Turn off plotting device
invisible(dev.off())
## Calculate how much time things took
end.time <- Sys.time()
time.taken <- difftime(end.time, start.time, units="mins")
cat("It took",time.taken,"minutes to complete this job\n")
}
