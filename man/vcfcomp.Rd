% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vcf-compare.R
\name{vcfcomp}
\alias{vcfcomp}
\title{Compare two VCF/BCF files reporting various statistics}
\usage{
vcfcomp(
  test,
  truth,
  formats = c("DS", "GT"),
  stats = "r2",
  by.sample = FALSE,
  by.variant = FALSE,
  flip = FALSE,
  names = NULL,
  bins = NULL,
  af = NULL,
  out = NULL,
  choose_random_start = FALSE,
  return_pse_sites = FALSE,
  ...
)
}
\arguments{
\item{test}{path to the comparision file (test), which can be a VCF/BCF file, vcftable object  or saved RDS file.}

\item{truth}{path to the baseline file (truth),which can be a VCF/BCF file, vcftable object  or saved RDS file.}

\item{formats}{character vector. the FORMAT tags to extract for the test and truth respectively.
default c("DS", "GT") extracts 'DS' of the target and 'GT' of the truth.}

\item{stats}{the statistics to be calculated. supports the following.
"r2": the Pearson correlation coefficient square.
"f1": the F1-score, good balance between sensitivity and precision.
"nrc": the Non-Reference Concordance rate
"pse": the Phasing Switch Error rate}

\item{by.sample}{logical. calculate sample-wise concordance, which can be stratified by MAF bin.}

\item{by.variant}{logical. calculate variant-wise concordance, which can be stratified by MAF bin.
If both by.sample and by.variant are FALSE, then do calculations for all samples and variants together in a bin.}

\item{flip}{logical. flip the ref and alt variants}

\item{names}{character vector. reset samples' names in the test VCF.}

\item{bins}{numeric vector. break statistics into allele frequency bins.}

\item{af}{file path with allele frequency or a RDS file with a saved object for af.
Format of the text file: a space-separated text file with five columns and a header named 'chr' 'pos' 'ref' 'alt' 'af'}

\item{out}{output prefix for saving objects into RDS file}

\item{choose_random_start}{choose random start for stats="pse"}

\item{return_pse_sites}{boolean. return phasing switch error sites}

\item{...}{options passed to \code{vcftable}}
}
\value{
a list of various statistics
}
\description{
Compare two VCF/BCF files reporting various statistics
}
\details{
\code{vcfcomp} implements various statisitcs to compare two VCF/BCF files,
e.g. report genotype concocrdance, correlation stratified by allele frequency.
}
\examples{
library('vcfppR')
test <- system.file("extdata", "imputed.gt.vcf.gz", package="vcfppR")
truth <- system.file("extdata", "imputed.gt.vcf.gz", package="vcfppR")
samples <- "HG00133,HG00143,HG00262"
res <- vcfcomp(test, truth, stats="f1", samples=samples, setid=TRUE)
str(res)
}
\author{
Zilong Li \email{zilong.dk@gmail.com}
}
