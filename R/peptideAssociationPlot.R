#' peptide-to-protein mapping
#' @param protein protein name.
#' @param batch batch No. in the Caprion experiment.
#' @param mm maximum number of mismatches.
#'
#' @export
#' @return a list containing mapping information.
#' @examples
#' \dontrun{
#'   batch <- "ZWK"
#'   load(paste0("~/Caprion/pilot/",batch,".rda"))
#'   PROC <- peptideMapping("PROC",mm=0)
#' }

peptideMapping <- function(protein,batch="ZWK",mm=5)
{
  Protein <- Modified.Peptide.Sequence <- Isotope.Group.ID <- Monoisotopic.mz <- Max.Isotope.Time.Centroid <- Charge <- NA
  accession <- subset(pQTLdata::caprion,grepl(protein,Protein))[["Accession"]]
  mapping <- subset(get(paste0("mapping_",batch)),grepl(protein,Protein))[1:6] %>%
             setNames(c("Isotope.Group.ID",
                        "Protein",
                        "Modified.Peptide.Sequence",
                        "Monoisotopic.mz",
                        "Max.Isotope.Time.Centroid",
                        "Charge"))
  mps <- mapping %>%
            dplyr::group_by(Modified.Peptide.Sequence) %>%
            dplyr::reframe(n_isotope=dplyr::n(),
                   isotope=paste(Isotope.Group.ID,collapse=";"),
                   mz=paste(Monoisotopic.mz,collapse=";"),
                   mtime=paste(Max.Isotope.Time.Centroid,collapse=";"),
                   charge=paste(Charge,collapse=";"))
  fasta_file_path <- paste0("https://rest.uniprot.org/uniprotkb/",accession,".fasta")
  fasta_sequences <- Biostrings::readAAStringSet(fasta_file_path, format = "fasta")
  sequence <- fasta_sequences[[1]]
  cat("Sequence:", toString(sequence), "\n")
  mp <- sapply(setdiff(mapping[["Modified.Peptide.Sequence"]],"-"),
        function(x) {
          y <- gsub("\\[\\d+\\.\\d+\\]", "X", x);
          mp <- Biostrings::matchPattern(y,sequence,with.indels=TRUE,max.mismatch=mm)
          as.data.frame(mp)
        })
  invisible(list(accession=accession,sequence=sequence,mapping=mapping,mps=mps,positions=t(mp)))
}

#' peptide association plot
#' @param protein protein name.
#' @param cistrans peptide signals and classification.
#' @param chrlen chromosome length (hg18, hg19, hg38).
#' @param disp distance from x-axis.
#'
#' @export
#' @return None.
#' @examples
#' \dontrun{
#'   par(mfrow=c(2,1))
#'   protein <- "PROC"
#'   suffix <- "_dr"
#'   input <- paste0("~/Caprion/analysis/METAL",suffix,"/gz/",protein,suffix,".txt.gz")
#'   annotation <- paste0("~/Caprion/analysis/METAL",suffix,"/vep/",protein,suffix,".txt")
#'   reference <- file.path(find.package("pQTLtools"),"turboman",
#'                                       "turboman_hg19_reference_data.rda")
#'   pvalue_sign <- 5e-8
#'   plot_title <- protein
#'   pQTLtools::turboman(input, annotation, reference, pvalue_sign, plot_title)
#'   cistrans <- read.csv(paste0("~/pQTLtools/tests","/",protein,".cis.vs.trans"))
#'   load("~/pQTLtools/tests/PROC.rda")
#'   peptideAssociationPlot(protein,cistrans)
#' }

peptideAssociationPlot <- function(protein,cistrans,chrlen=gap::hg19,disp=85)
{
  SNPChrom <- SNPPos <- isotope <- NA
  par(mar=c(16,3,1,1))
  mapping <- get(protein)
  g2d <-  gap::grid2d(chrlen,plot=FALSE)
  n <- with(g2d, n-1)
  CM <- with(g2d, CM)
  dseq <- CM[n+1]/length(mapping$sequence)
  positions <- with(mapping,positions) %>%
               data.frame
  positions <- data.frame(positions,Modified.Peptide.Sequence=rownames(positions),ID=(1:nrow(positions))/2)
  cistrans <- mutate(cistrans,pos=g2d$CM[SNPChrom]+SNPPos) %>%
              left_join(mapping$mapping,by=c('isotope'='Isotope.Group.ID')) %>%
              left_join(positions)
  chr <- cistrans[["SNPChrom"]]
  chr[chr == "X"] <- 23
  chr[chr == "Y"] <- 24
  pos <- CM[chr] + cistrans[["SNPPos"]]
  xy <- xy.coords(c(0, CM), seq(0, max(cistrans[["log10p"]]),length.out=n+3))
  par(xaxt = "n", yaxt = "n")
  plot(xy$x, xy$y, type = "n", ann = FALSE, axes = FALSE)
  par(xaxt = "s", yaxt = "s", xpd = TRUE)
  axis(2, pos = 2, cex = 1.2)
  title(xlab="", ylab = "-log10(P)", line = 1)
  all_colors <- colors()
  all_numbers <- 1:length(all_colors)
  nonwhite <- !grepl("white",all_colors)
  for(iso in unique(cistrans[["isotope"]]))
  {
    d <- filter(cistrans,isotope==iso)
    c <- d[["SNPChrom"]]
    p <- CM[c]+d[["SNPPos"]]
    points(p, d[["log10p"]], cex = 0.8, col = ifelse(d[["cis"]], "red", "blue"), pch = 19)
    lines(c(as.integer(d[["start"]]), as.integer(d[["end"]]))*dseq, -c(d[["ID"]]+disp, d[["ID"]]+disp),
            col=all_numbers[nonwhite][d[["ID"]]], lwd = 4)
    segments((as.integer(d[["start"]])+as.integer(d[["end"]]))*dseq/2,-(d[["ID"]]+disp), p, d[["log10p"]]*1,
            col=all_numbers[nonwhite][d[["ID"]]])
  }
  for (x in 1:n) {
       segments(CM[x], 0, CM[x], max(cistrans[["log10p"]]), col = "black")
       text(ifelse(x == 1, CM[x+1]/2, (CM[x+1] + CM[x])/2),
            0, pos = 1, offset = 0.5, gap::xy(x), cex = 0.8)

  }
  segments(0, 0, CM[n+1], 0)
}
