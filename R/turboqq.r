#' QQ plots
#'
#' @param input_data_path Path of the input association data.
#' @param output_data_rootname Root name of the plot output file.
#' @param plot_title Plot title to be displayed on top of the plot.
#' @export
#' @return No direct return value. The script generates QQ plots as output.
#'
#' @details
#' The method uses the kth order statistic from a sample of n i.i.d. U(0,1) statistics has a Beta(k,n+1-k) distribution as in
#' \insertCite{quesenberry80;textual}{pQTLtools} and Coded by  Weale M, Price T.
#' \href{https://sites.google.com/site/mikeweale/software}{https://sites.google.com/site/mikeweale/software}
#'
#' **Input association data path / input_data_path:**
#'
#' Define path of the input association data. The input data needs to be a file that has:
#' 1. Spaces as field separators
#' 2. One header line
#' 3. Option I. (no extreme p-values present): 3 columns, being chromosome, position, pvalue in order, column names are not important.
#'    Option II. (extreme p-values present): 5 columns, being chromosome, position, pvalue, beta, se in order, column names are not important.
#'
#' @examples
#' \dontrun{
#' require(gap.datasets)
#' test <- mhtdata[c("chr","pos","p")]
#' write.table(test,file="test.txt",row.names=FALSE,quote=FALSE)
#' input_data_path <- "test.txt"
#' output_data_rootname <- "test_qq"
#' plot_title <- "gap.datasets example"
#' turboqq(input_data_path, output_data_rootname, plot_title)
#' }
#'
#' @references
#' \insertAllCited{}
#'
#' @author Bram Prins

turboqq <- function(input_data_path, output_data_rootname, plot_title) {
  
  # Verbose printing status bars
  fat_status_bar <- "============================================================================================================"
  skinny_status_bar <- "------------------------------------------------------------------------------------------------------------"
  
  cat("\n", fat_status_bar, "\n 1. Read in arguments from the command line\n", fat_status_bar, "\n\n")
  cat(paste0("  Data file path : ", input_data_path), "\n")
  
  # Assign variable classes
  input_data_path <- as.character(input_data_path)
  output_data_rootname <- as.character(output_data_rootname)
  plot_title <- as.character(plot_title)
  
  # Reading in association plot data with scan
  cat("\n", fat_status_bar, "\n 2. Reading in association plot data with scan\n", fat_status_bar, "\n\n")
  
  initial_data_dims <- dim(as.data.frame(read.table(input_data_path, header=TRUE, stringsAsFactors=FALSE, nrows=10)))[2]
  
  if (initial_data_dims == 3) {
    initial_data <- data.frame(scan(input_data_path,
                                    what = list(chromosome = 0, position = 0, pvalue = 0),
                                    skip = 1,
                                    sep = " ",
                                    quiet = TRUE))
    initial_data_contains_beta_se <- FALSE
  } else if (initial_data_dims == 5) {
    initial_data <- data.frame(scan(input_data_path,
                                    what = list(chromosome = 0, position = 0, pvalue = 0, beta = 0, se = 0),
                                    skip = 1,
                                    sep = " ",
                                    quiet = TRUE))
    initial_data_contains_beta_se <- TRUE
  } else {
    stop("Input data does not have expected dimensions")
  }
  
  cat("\n", fat_status_bar, "\n 3. Calculate log P values\n", fat_status_bar, "\n\n")
  
  # Check if p-values are already logged
  if (length(which(initial_data$pvalue > 1)) > 0) {
    initial_data$log_pvalue <- initial_data$pvalue
    # Get only the complete data
    initial_data <- initial_data[complete.cases(initial_data), ]
    # Remove the original pvalues
    initial_data$pvalue <- NULL
  } else {
    # Calculate the -log10 p-value for the input data
    initial_data$log_pvalue <- -log10(initial_data$pvalue)
    # If beta/SE are provided, and pvalues are missing (because they are extreme), log10 P recalculate from beta/SE
    missing_pvalues_index <- which((is.na(initial_data$log_pvalue) | initial_data$log_pvalue == 0))
    
    if (initial_data_contains_beta_se & (length(missing_pvalues_index) > 0)) {
      # Calculate expected p-values for missing data
      missing_pvalues <- (-log(2, base = 10) - pnorm(-abs(initial_data$beta[missing_pvalues_index]/initial_data$se[missing_pvalues_index]), log = T) / log(10))
      # Only replace if indeed they were below the smallest non-zero normalized floating-point number
      initial_data[missing_pvalues_index, c("log_pvalue")] <- ifelse(missing_pvalues > -log10(.Machine$double.xmin), missing_pvalues, NA)
    }
    # Get only the complete data
    initial_data <- initial_data[complete.cases(initial_data), ]
    # Remove the original pvalues
    initial_data$pvalue <- NULL
  }
  
  plot_resolution <- 1800
  
  # For calculating concentration bands
  alpha <- 0.05
  df <- 1
  one.sided <- FALSE
  
  # Calculate plotting parameters
  cat("\n", fat_status_bar, "\n 4. Calculate plotting points\n", fat_status_bar, "\n\n")
  
  # Sort pvals and calculate expected pvals
  log_obspval_sorted <- sort(initial_data$log_pvalue, decreasing = TRUE)
  exppval <- c(1:length(log_obspval_sorted))
  log_exppval <- -(log10((exppval - 0.5) / length(exppval)))
  cat(unique(sort(round(log_exppval, 0))), "\n")
  
  # Get the maxima
  ymax <- max(na.omit(log_obspval_sorted))
  xmax <- max(na.omit(log_exppval))
  
  # Calculate the lambda
  lambda <- median(qchisq(na.omit((-log_obspval_sorted * log(10))), df = 1, lower.tail = FALSE, log.p = TRUE)) / qchisq(0.5, 1)
  if (!is.null(initial_data$beta) & !is.null(initial_data$se))
    lambda <- median((initial_data$beta / initial_data$se)^2, na.rm = TRUE) / qchisq(0.5, 1)
  
  # Set vertical resolution, allowing a fixed number of points to be plotted vertically
  vertical_resolution <- plot_resolution
  horizontal_resolution <- plot_resolution
  
  # Scale the resolution for the p-values
  obs_log_pvalue_break_size <- ymax / vertical_resolution
  exp_log_pvalue_break_size <- xmax / horizontal_resolution
  
  # Create a vector from 0 to the vertical resolution, which we will use to bin pvalues
  obs_log_pvalue_scaling_vector <- seq(0, vertical_resolution, by = obs_log_pvalue_break_size)
  exp_log_pvalue_scaling_vector <- seq(0, horizontal_resolution, by = exp_log_pvalue_break_size)
  
  # Bin the pvals and get only unique pairs
  obs_log_pvalue_binned <- .bincode(log_obspval_sorted, obs_log_pvalue_scaling_vector, right = TRUE, include.lowest = FALSE) * obs_log_pvalue_break_size
  exp_log_pvalue_binned <- .bincode(log_exppval, exp_log_pvalue_scaling_vector, right = TRUE, include.lowest = FALSE) * exp_log_pvalue_break_size
  plot_data_binned <- as.data.frame(cbind(exp_log_pvalue_binned, obs_log_pvalue_binned))
  plot_data_reduced <- unique(plot_data_binned)
  
  # Truncate the maxima and use these to make nice ticks
  pretty_ymax = trunc(ymax + 1)
  pretty_bigymax = trunc(ymax)
  pretty_xmax = trunc(xmax + 1)
  
  y4Ly <- switch(
    TRUE,
    pretty_ymax = c(0:pretty_ymax),
    (pretty_ymax <= 15) ~ c(seq(0, pretty_bigymax + 10, 10)),
    (pretty_ymax <= 100) ~ c(seq(0, pretty_bigymax + 20, 20)),
    (pretty_ymax <= 200) ~ c(seq(0, pretty_bigymax + 30, 30)),
    (pretty_ymax <= 300) ~ c(seq(0, pretty_bigymax + 40, 40)),
    (pretty_ymax <= 400) ~ c(seq(0, pretty_bigymax + 50, 50)),
    (pretty_ymax <= 500) ~ c(seq(0, pretty_bigymax + 60, 60)),
    (pretty_ymax <= 600) ~ c(seq(0, pretty_bigymax + 70, 70)),
    (pretty_ymax <= 700) ~ c(seq(0, pretty_bigymax + 80, 80)),
    (pretty_ymax <= 800) ~ c(seq(0, pretty_bigymax + 90, 90)),
    (pretty_ymax <= 900) ~ c(seq(0, pretty_bigymax + 100, 100)),
    TRUE ~ c(seq(0, pretty_bigymax + 200, 200))
  )
  
  x4Lx = c(0:pretty_xmax)
  xnums = (x4Lx)
  ynums = (y4Ly)
  Lx <- parse(text = paste(x4Lx, sep = ""))
  Ly <- parse(text = paste(y4Ly, sep = ""))
  
  # Use the maxima for making x and y limits for plotting function
  x.lim = xmax
  y.lim = max(y4Ly)
  
  # Start the actual plotting
  cat("\n", fat_status_bar, "\n 5. Start the actual plotting\n", fat_status_bar, "\n\n")
  
  # Start device
  png(paste0(output_data_rootname, ".png"), height = plot_resolution, width = plot_resolution, pointsize = 12, res = 450)
  par(mar = c(4, 4, 3, 1))
  
  # Create empty plot
  # Just plots the outside box
  plot(0, main = plot_title, xlab = "Expected p-value", ylab = "Observed p-value", type = "n", xlim = c(0, x.lim), ylim = c(0, y.lim), xaxt = 'n', yaxt = 'n')
  # Draw axes
  axis(1, at = xnums, labels = unique(sort(round(log_exppval, 0)))[1:length(xnums)], las = 2)
  axis(2, at = ynums, labels = Ly, las = 2)
  
  # Define and draw the concentration bands
  # Note that conc band won't draw if x has too many datapoints
  n <- length(log_exppval)
  frac = 1
  starti = floor((n - 1) * (1 - frac)) + 1
  lena = n - starti + 1
  a2 = (1:lena)
  b <- n + 1 - a2
  
  cat("\n", fat_status_bar, "\n  a. Calculate and plot concentration band points\n", fat_status_bar, "\n\n")
  
  # Define and draw the concentration bands
  # Note that conc band won't draw if x has too many datapoints
  if (one.sided == FALSE) {
    upper <- -log10(qbeta(1 - alpha/2, a2, b)) # Exp. upper CL for 'a'th U(0,1) order statistic (becomes 'lower')
    lower <- -log10(qbeta(alpha/2, a2, b)) # Exp. lower CL for 'a'th U(0,1) order statistic (becomes 'upper')
  } else {
    upper <- rep(1, length(log_exppval)) # Exp. upper CL for 'a'th U(0,1) order statistic (becomes 'lower')
    lower <- qbeta(alpha, a2, b) # Exp. lower CL for 'a'th U(0,1) order statistic (becomes 'upper')
  }
  
  polygon(c(log_exppval, rev(log_exppval)), c(upper, rev(log_exppval)), col = "grey", border = NA) # 'lower' band
  polygon(c(log_exppval, rev(log_exppval)), c(lower, rev(log_exppval)), col = "grey", border = NA) # 'upper' band
  
  # Draw the diagonal
  lines(c(0, xmax), c(0, xmax), col = "blue", lty = 2, lwd = 0.5)
  
  # Print the lambda
  legend("topleft", legend = substitute(paste(lambda[GC], "=", x), list(x = formatC(round(lambda, 3), 3, format = "f"))),
         text.col = ifelse(lambda > 1.1, "red", "black"), bty = "n")
  
  cat("\n", skinny_status_bar, "\n  b. Plot points P values\n", skinny_status_bar, "\n\n")
  
  # Plot points
  points(plot_data_reduced[, 1], plot_data_reduced[, 2], pch = 19, cex = 0.5, col = "dodgerblue4")
  
  # Turn off the device
  dev.off()
}
