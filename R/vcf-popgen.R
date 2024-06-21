#' @title
#' count the heterozygous sites per sample in the VCF/BCF
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
#' @param fun which popgen function to run. available functions are
#'            "heterozygosity".
#'
#' @return \code{vcfpopgen} a list containing the following components:
#'\describe{
#'\item{samples}{: character vector; \cr
#'                 the samples ids in the VCF file after subsetting
#'}
#' 
#'\item{hets}{: integer vector; \cr
#'              the counts of heterozygous sites of each sample in the same order as \code{samples}
#'}
#'
#'}
#' @author Zilong Li \email{zilong.dk@gmail.com}
#'
#' @examples
#' library('vcfppR')
#' vcffile <- system.file("extdata", "raw.gt.vcf.gz", package="vcfppR")
#' res <- vcfpopgen(vcffile)
#' str(res)
#' @export
vcfpopgen <- function(vcffile,
                      region = "",
                      samples = "-",
                      pass = FALSE,
                      qual = 0,
                      fun = "heterozygosity") {
  if(!file.exists(vcffile))
    stop("file doesn't exist")
  return(heterozygosity(vcffile, region, samples, pass, qual))
}
