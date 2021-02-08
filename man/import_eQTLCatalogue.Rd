\name{import_eQTLCatalogue}
\alias{import_eQTLCatalogue}
\title{Import eQTL Catalogue}
\usage{import_eQTLCatalogue(ftp_path, region, selected_gene_id, column_names, verbose = TRUE)}
\description{
An adopted function which imports eQTL Catalogue.
}
\arguments{
\item{ftp_path}{URL}
\item{region}{chr:start-end}
\item{selected_gene_id}{An Ensembl gene ID}
\item{column_names}{Column names of the dataset}
\item{verbose}{Extra information}
}
\details{
This function is based on the eQTL-Catalogue-resources.
}

\value{
A summary statistic object
}

\section{References}{
Kerimov, N., et al. (2020). "eQTL Catalogue: a compendium of uniformly processed human gene expression and splicing QTLs." bioRxiv: 2020.2001.2029.924266.
	https://www.ebi.ac.uk/eqtl/.
}

\examples{
\dontrun{
library(pQTLtools)
invisible(lapply(c("dplyr", "ggplot2", "readr", "coloc", "GenomicRanges","seqminer"), require, character.only = TRUE))
tfp <- file.path(path.package("pQTLtools"),"eQTL-Catalogue","tabix_ftp_paths.tsv")
tabix_paths <- read.delim(tfp, sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>% dplyr::as_tibble()
tfpi <- file.path(path.package("pQTLtools"),"eQTL-Catalogue","tabix_ftp_paths_imported.tsv")
imported_tabix_paths <- read.delim(tfpi, sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>% dplyr::as_tibble()

# MPV association at the ARHGEF3 locus
region <- "3:56615721-57015721"
ensGene <- "ENSG00000163947"
platelet_df <- dplyr::filter(tabix_paths, study == "CEDAR", tissue_label == "platelet")
hdr <- file.path(path.package("pQTLtools"),"eQTL-Catalogue","column_names.CEDAR")
column_names <- names(read.delim(hdr))
summary_stats <- import_eQTLCatalogue(platelet_df$ftp_path, region, selected_gene_id = ensGene, column_names)
summary_stats
ggplot(summary_stats, aes(x = position, y = -log(pvalue, 10))) + geom_point()
# gwasvcf::set_bcftools(path=file.path(HPC_WORK,"bin","bcftools"))
# GWAS sumstat from the same region
# manually download and parse with gwasvcf
# wget https://gwas.mrcieu.ac.uk/files/ebi-a-GCST004599/ebi-a-GCST004599.vcf.gz
# wget https://gwas.mrcieu.ac.uk/files/ebi-a-GCST004599/ebi-a-GCST004599.vcf.gz.tbi
# gwas_stats <- gwasvcf::query_gwas("ebi-a-GCST004599.vcf.gz", chrompos = "3:56649749-57049749")
# gwas_stats <- gwasvcf::vcf_to_granges(gwas_stats) %>% keepSeqlevels("3") %>% renameSeqlevels("chr3")
# via import_OpenGWAS
opengwas_id <- "ebi-a-GCST004599"
region <- "3:56649749-57049749"
gwas_stats <- import_OpenGWAS(opengwas_id,region) %>% keepSeqlevels("3") %>% renameSeqlevels("chr3")
f <- file.path(path.package("pQTLtools"),"eQTL-Catalogue","hg19ToHg38.over.chain")
chain <- rtracklayer::import.chain(f)
gwas_stats_hg38 <- rtracklayer::liftOver(gwas_stats, chain) %>%
  unlist() %>%
  renameSeqlevels("3") %>%
  dplyr::as_tibble() %>%
  dplyr::transmute(chromosome = seqnames, position = start, AF, ES, SE, LP, SS) %>%
  dplyr::mutate(id = paste(chromosome, position, sep = ":")) %>%
  dplyr::mutate(MAF = pmin(AF, 1-AF)) %>%
  dplyr::group_by(id) %>%
  dplyr::mutate(row_count = n()) %>%
  dplyr::ungroup() %>%
  dplyr::filter(row_count == 1)
ggplot(gwas_stats_hg38, aes(x = position, y = LP)) + geom_point()
# Colocalisation
res <- run_coloc(summary_stats, gwas_stats_hg38)

# a. all other eQTL datasets
microarray_df <- dplyr::filter(tabix_paths, quant_method == "microarray") %>% dplyr::mutate(qtl_id = paste(study, qtl_group, sep = "_"))
ftp_path_list <- setNames(as.list(microarray_df$ftp_path), microarray_df$qtl_id[1])
hdr <- file.path(path.package("pQTLtools"),"eQTL-Catalogue","column_names.CEDAR")
column_names <- names(read.delim(hdr))
summary_list <- purrr::map(ftp_path_list, ~import_eQTLCatalogue(., region, selected_gene_id = ensGene, column_names))
coloc_df_microarray <- purrr::map_df(summary_list, ~run_coloc(., gwas_stats_hg38), .id = "qtl_id")

# b. Uniformly processed RNA-seq datasets
rnaseq_df <- dplyr::filter(tabix_paths, quant_method == "ge") %>% dplyr::mutate(qtl_id = paste(study, qtl_group, sep = "_"))
ftp_path_list <- setNames(as.list(rnaseq_df$ftp_path), rnaseq_df$qtl_id)
hdr <- file.path(path.package("pQTLtools"),"eQTL-Catalogue","column_names.Alasoo")
column_names <- names(read.delim(hdr))
safe_import <- purrr::safely(import_eQTLCatalogue)
summary_list <- purrr::map(ftp_path_list, ~safe_import(., region, selected_gene_id = ensGene, column_names))
result_list <- purrr::map(summary_list, ~.$result)
result_list <- result_list[!unlist(purrr::map(result_list, is.null))]
coloc_df_rnaseq <- purrr::map_df(result_list, ~run_coloc(., gwas_stats_hg38), .id = "qtl_id")

# c. GTEx_v8 imported eQTL datasets
rnaseq_df <- dplyr::filter(imported_tabix_paths, quant_method == "ge") %>% dplyr::mutate(qtl_id = paste(study, qtl_group, sep = "_"))
ftp_path_list <- setNames(as.list(rnaseq_df$ftp_path), rnaseq_df$qtl_id)
hdr <- file.path(path.package("pQTLtools"),"eQTL-Catalogue","column_names.GTEx")
column_names <- names(read.delim(hdr))
safe_import <- purrr::safely(import_eQTLCatalogue)
summary_list <- purrr::map(ftp_path_list, ~safe_import(., region, selected_gene_id = ensGene, column_names))
result_list <- purrr::map(summary_list, ~.$result)
result_list <- result_list[!unlist(purrr::map(result_list, is.null))]
result_filtered <- purrr::map(result_list, ~dplyr::filter(., !is.na(se)))
coloc_df_imported <- purrr::map_df(result_filtered, ~run_coloc(., gwas_stats_hg38), .id = "qtl_id")

coloc_df = dplyr::bind_rows(coloc_df_microarray, coloc_df_rnaseq, coloc_df_imported)
dplyr::arrange(coloc_df, -PP.H4.abf)
ggplot(coloc_df, aes(x = PP.H4.abf)) + geom_histogram()
}
}
\note{Adapted function}
\keyword{utilities}
