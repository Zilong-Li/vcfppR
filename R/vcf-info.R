#' @title
#' read a INFO tag in the VCF/BCF into R data structure
#'
#' @details
#' \code{vcftable} uses the C++ API of vcfpp, which is a wrapper of htslib, to read VCF/BCF files.
#' Thus, it has the full functionalities of htslib, such as restrict to specific variant types,
#' samples and regions. For the memory efficiency reason, the \code{vcftable} is designed
#' to parse only one tag at a time in the FORMAT column of the VCF. In default, only the matrix of genotypes,
#' i.e. "GT" tag, are returned by \code{vcftable}, but there are many other tags supported by the \code{format} option.
#'
#' @param vcffile path to the VCF/BCF file
#'
#' @param tag the INFO tag to extract.
#'
#' @param region region to subset in bcftools-like style: "chr1", "chr1:1-10000000"
#'
#' @param vartype restrict to specific type of variants. supports "snps","indels", "sv", "multisnps","multiallelics"
#' 
#' @param ids  character vector. restrict to sites with ID in the given vector. default NULL won't filter any sites.
#'
#' @param qual numeric. restrict to variants with QUAL > qual.
#'
#' @param pass logical. restrict to variants with FILTER = "PASS".
#'
#' @param setid logical. reset ID column as CHR_POS_REF_ALT.
#'
#' @return Return a list containing the following components:
#'\describe{
#'\item{chr}{: character vector; \cr
#'             the CHR column in the VCF file
#'}
#'
#'\item{pos}{: character vector; \cr
#'             the POS column in the VCF file
#'}
#'
#'\item{id}{: character vector; \cr
#'            the ID column in the VCF file
#'}
#'
#'\item{ref}{: character vector; \cr
#'             the REF column in the VCF file
#'}
#'
#'\item{alt}{: character vector; \cr
#'             the ALT column in the VCF file
#'}
#'
#'\item{qual}{: character vector; \cr
#'             the QUAL column in the VCF file
#'}
#'
#'\item{filter}{: character vector; \cr
#'                the FILTER column in the VCF file
#'}
#'
#'\item{tag}{: vector of either integer, numberic or character values depending on the tag to extract; \cr
#'             a specifiy tag in the INFO column to be extracted
#'}
#'}
#' @author Zilong Li \email{zilong.dk@gmail.com}
#'
#' @examples
#' library('vcfppR')
#' vcffile <- system.file("extdata", "raw.gt.vcf.gz", package="vcfppR")
#' res <- vcfinfo(vcffile, "AF", region = "chr21:1-5050000", vartype = "snps")
#' str(res)
#
#' @export
vcfinfo <- function(vcffile,
                    tag,
                    region = "",
                    vartype = "all",
                    ids = NULL,
                    qual = 0,
                    pass = FALSE,
                    setid = FALSE) {
  snps <- FALSE
  indels <- FALSE
  svs <- FALSE
  multiallelics <- FALSE
  multisnps <- FALSE
  if(vartype == "snps") snps <- TRUE
  else if(vartype == "indels") indels <- TRUE
  else if(vartype == "sv") svs <- TRUE
  else if(vartype == "multisnps") multisnps <- TRUE
  else if(vartype == "multiallelics") multiallelics <- TRUE
  else if(vartype != "all") stop("Invaild variant type!")
  if(is.null(ids)) ids <- c("")
  res <- tableInfo(vcffile, tag, region, ids, qual, pass, snps, indels, multiallelics, multisnps, svs)
  return(res)
}
