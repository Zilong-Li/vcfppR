#' @title
#' read VCF/BCF contents into R data structure
#'
#' @description
#' The swiss army knife for reading VCF/BCF into R data types rapidly and easily.
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
#' @param region region to subset in bcftools-like style: "chr1", "chr1:1-10000000"
#'
#' @param samples samples to subset in bcftools-like style.
#'                comma separated list of samples to include (or exclude with "^" prefix).
#'                e.g. "id01,id02", "^id01,id02".
#'
#' @param vartype restrict to specific type of variants. supports "snps","indels", "multisnps","multiallelics"
#' @param qual logical. restrict to variants with QUAL > qual.
#'
#' @param pass logical. restrict to variants with FILTER = "PASS".
#'
#' @param info logical. drop INFO column in the returned list.
#'
#' @param format the FORMAT tag to extract.
#'
#' @param collapse logical. It acts on the FORMAT. If the FORMAT is GT, the dim of raw genotypes matrix of diploid is (M, 2 * N),
#'                 where M is #markers and N is #samples. default TRUE will collapse the genotypes for each sample such that the matrix is (M, N).
#'                 set this to false if one wants to maintain the phasing order, e.g. "1|0" is parsed as c(1, 0) with collapse=FALSE.
#'                 If the FORMAT is not GT, then it will turn a list of the format vector into a matrix with collapse=TRUE. However, this raises issues
#'                 when the variant is mutliallelic because it has more vaules for this row.
#'
#' @return \code{vcftable} a list containing the following components:
#'\describe{
#'\item{samples}{: character vector; \cr
#'                 the samples ids in the VCF file after subsetting
#'}
#'
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
#'\item{info}{: character vector; \cr
#'              the INFO column in the VCF file
#'}
#'
#'\item{format}{: matrix of either integer of numberic values depending on the tag to extract; \cr
#'                a specifiy tag in the FORMAT column to be extracted
#'}
#'}
#' @author Zilong Li \email{zilong.dk@gmail.com}
#'
#' @examples
#' library('vcfppR')
#' vcffile <- system.file("extdata", "raw.gt.vcf.gz", package="vcfppR")
#' res <- vcftable(vcffile, "chr21:1-5050000", vartype = "snps")
#' str(res)
#' @export
vcftable <- function(vcffile, region = "", samples = "-", vartype = "all", qual = 0,
                     pass = FALSE, info = TRUE, format = "GT", collapse = TRUE) {
  snps <- FALSE
  indels <- FALSE
  multiallelics <- FALSE
  multisnps <- FALSE
  if(vartype == "snps") snps <- TRUE
  else if(vartype == "indels") indels <- TRUE
  else if(vartype == "multisnps") multisnps <- TRUE
  else if(vartype == "multiallelics") multiallelics <- TRUE
  else if(vartype != "all") stop("Invaild variant type!")
  res <- NULL
  if(format == "GT") {
    res <- tableGT(vcffile, region, samples, qual, pass, info, snps, indels, multiallelics, multisnps)
    res[[10]] <- do.call("rbind", res[[10]])
    n <- ncol(res$gt)
    ploidy <- n / length(res$samples)
    if(ploidy == 2 && collapse) {
      res$gt <- res$gt[, seq(1, n, 2)] + res$gt[, seq(2, n, 2)]
      res$gt[res$gt < 0] <- NA
    } else {
      res$gt[res$gt < 0] <- NA
    }
  } else {
    res <- tableOther(format, vcffile, region, samples, qual, pass, info, snps, indels, multiallelics, multisnps)
    if(collapse) res[[10]] <- do.call("rbind", res[[10]])
  }
  return(res)
}
