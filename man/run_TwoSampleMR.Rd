\name{run_TwoSampleMR}
\alias{run_TwoSampleMR}
\title{A generic wrapper for TwoSampleMR analysis}
\usage{run_TwoSampleMR(exposure.args=format_data.args(),outcome.args=extract_outcome_data.args(),prefix,...))}
\description{
Given a compressed or uncompressed file containing exposure data, this function conducts a two-sample MR analysis using MR-Base.
}
\arguments{
\item{exposure.args}{a list of options appropriate for TwoSampleMR's format_data()}
\item{outcome.args}{a list of options appropriate for TwoSampleMR's extract_outcome_data()}
\item{prefix}{a prefix for output files}
\item{...}{other options to be added}
\item{file}{file}
\item{type}{exposure or outcome}
\item{phenotype_col}{protein column}
\item{header}{TRUE or FALSE}
\item{snp_col}{rsid}
\item{effect_allele_col}{effect allele}
\item{other_allele_col}{other allele}
\item{eaf_col}{EAF}
\item{beta_col}{beta}
\item{se_col}{se}
\item{pval_col}{pval}
\item{log_pval}{FALSE or TRUE}
\item{snp}{snpname}
\item{outcomes}{outcome}
\item{proxies}{TRUE or FALSE}
\item{rsq}{r^2}
\item{align_alleles}{1}
\item{palindromes}{1}
\item{maf_threshold}{0.3}
\item{access_token}{token}
\item{splitsize}{10000}
\item{proxy_splitsize}{500}
}
\details{
As TwoSampleMR includes seemingly daunting tasks, this function intends to simplify various steps in a two-sample MR. It is
particularly useful when a large numbher of MRs are necessary, e.g., multiple proteins and their cis/trans regions need to be examined,
in which case prefix could direct the output to various directories.
}

\value{
No value is returned but several files.
}

\section{References}{
Dimou NL, Tsilidis KK. A Primer in Mendelian Randomization Methodology with a Focus on Utilizing Published Summary Association Data. In
Evangelos Evangelou (ed.), Genetic Epidemiology: Methods and Protocols, Methods in Molecular Biology, vol. 1793,
https://doi.org/10.1007/978-1-4939-7868-7_13, Springer Science+Business Media, LLC, part of Springer Nature 2018
}

\examples{
\dontrun{
To be added.
}
}
\note{adapted from script by Dimou NL, Tsilidis KK}
\keyword{utilities}
