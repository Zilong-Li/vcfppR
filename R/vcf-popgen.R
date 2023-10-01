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
#' @return \code{popgen.heterozygosity} a list containing the following components:
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
#' @export
popgen.heterozygosity <- function(vcffile, region = "", samples = "-", pass = FALSE, qual = 0) {
    return(heterozygosity(vcffile, region, samples, pass, qual))
}
