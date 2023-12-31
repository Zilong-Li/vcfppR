% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vcf-popgen.R
\name{vcfpopgen}
\alias{vcfpopgen}
\title{count the heterozygous sites per sample in the VCF/BCF}
\usage{
vcfpopgen(
  vcffile,
  region = "",
  samples = "-",
  pass = FALSE,
  qual = 0,
  fun = "heterozygosity"
)
}
\arguments{
\item{vcffile}{path to the VCF/BCF file}

\item{region}{region to subset like bcftools}

\item{samples}{samples to subset like bcftools}

\item{pass}{restrict to variants with FILTER==PASS}

\item{qual}{restrict to variants with QUAL > qual.}

\item{fun}{which popgen function to run. available functions are
"heterozygosity".}
}
\value{
\code{vcfpopgen} a list containing the following components:
\describe{
\item{samples}{: character vector; \cr
the samples ids in the VCF file after subsetting
}

\item{hets}{: integer vector; \cr
the counts of heterozygous sites of each sample in the same order as \code{samples}
}

}
}
\description{
count the heterozygous sites per sample in the VCF/BCF
}
\examples{
library('vcfppR')
vcffile <- system.file("extdata", "raw.gt.vcf.gz", package="vcfppR")
res <- vcfpopgen(vcffile)
str(res)
}
\author{
Zilong Li \email{zilong.dk@gmail.com}
}
