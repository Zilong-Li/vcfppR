#' @title
#' summarize the various variant types at both variant level and sample level.
#'
#' @details
#' bcftools view -s "id01,id02" input.bcf.gz chr1:100000-20000
#'
#' @param vcffile path to the VCF/BCF file
#'
#' @param region region to subset like bcftools
#'
#' @param samples samples to subset like bcftools
#'
#' @param pass restrict to variants with FILTER==PASS
#'
#' @param qual restrict to variants with QUAL > qual.
#'
#' @param svtype summarize the variants with SVTYPE
#'
#' @return \code{vcfsummary} a list containing the following components:
#'\describe{
#'\item{summary}{: named integer vector; \cr
#'                 summarize the counts of each variant type
#'}
#'
#'\item{samples}{: character vector; \cr
#'                 the samples ids in the VCF file after subsetting
#'}
#' 
#'\item{vartype}{: integer vector; \cr
#'                 the counts of the variant type at sample level in the same order as \code{samples}
#'}
#' 
#'}
#' @author Zilong Li \email{zilong.dk@gmail.com}
#'
#' @examples
#' library('vcfppR')
#' svfile <- system.file("extdata", "sv.vcf.gz", package="vcfppR")
#' res <- vcfsummary(svfile, region = "chr21:1-10000000", svtype = TRUE)
#' str(res)
#' @export
vcfsummary <- function(vcffile,
                       region = "",
                       samples = "-",
                       pass = FALSE,
                       qual = 0,
                       svtype = FALSE) {
  if(!file.exists(vcffile))
    stop("file doesn't exist")
  if(svtype) {
    ret <- summarySVs(vcffile, region, samples, pass, qual)
  } else {
    ret <- summaryVariants(vcffile, region, samples, pass, qual)
  }
  class(ret) <- "vcfsummary"
  return(ret)
}
