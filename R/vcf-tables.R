#' @title
#' read VCF/BCF contents into R data structure
#'
#' @description
#' If you want to read VCF/BCF into R data types fast and simplely, this is the swiss army knife
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
#' @param format the FORMAT tag to extract, valid values are
#'               "GT","GP", "DP","DS","GL","PL","GQ","HQ","MQ","PQ"
#'
#' @param ploidy the ploidy of the organism if it is not diploidy
#'
#' @return \code{vcftable} a list containing the following components:
#'\describe{
#'\item{samples}{  character vector; \cr
#'                 the samples ids in the VCF file after subsetting
#'}
#'
#'\item{chr}{  character vector; \cr
#'             the CHR column in the VCF file
#'}
#'
#'\item{pos}{  character vector; \cr
#'             the POS column in the VCF file
#'}
#'
#'\item{id}{  character vector; \cr
#'            the ID column in the VCF file
#'}
#'
#'\item{ref}{  character vector; \cr
#'             the REF column in the VCF file
#'}
#'
#'\item{alt}{  character vector; \cr
#'             the ALT column in the VCF file
#'}
#'
#'\item{qual}{  character vector; \cr
#'             the QUAL column in the VCF file
#'}
#'
#'\item{filter}{  character vector; \cr
#'                the FILTER column in the VCF file
#'}
#'
#'\item{info}{  character vector; \cr
#'              the INFO column in the VCF file
#'}
#'
#'\item{format}{  matrix of either integer of numberic values depending on the tag to extract; \cr
#'                a specifiy tag in the FORMAT column to be extracted
#'}
#'}
#' @author Zilong Li \email{zilong.dk@gmail.com}
#'
#' @examples
#' library('vcfppR')
#' vcffile <- system.file("extdata", "test-GL.vcf.gz", package="vcfppR")
#' (res <- vcftable(vcffile,"chr20"))
#' vcffile <- system.file("extdata", "test-PL.vcf.gz", package="vcfppR")
#' (res <- vcftable(vcffile,"chr20", format = "PL"))
#' str(res)
#' @export
vcftable <- function(vcffile, region = "", samples = "-", format = "GT", ploidy = 2) {
  res <- switch(format,
                GT = tableGT(vcffile, region, samples),
                GP = tableGP(vcffile, region, samples),
                GQ = tableGQ(vcffile, region, samples),
                DS = tableDS(vcffile, region, samples),
                DP = tableDP(vcffile, region, samples),
                PL = tablePL(vcffile, region, samples),
                GL = tableGL(vcffile, region, samples),
                HQ = tableHQ(vcffile, region, samples),
                MQ = tableMQ(vcffile, region, samples),
                PQ = tablePQ(vcffile, region, samples),
                stop("Invaild tag in FORAMT column"))
  res[[10]] <- do.call("rbind", res[[10]])
  if(ploidy == 2 && format == "GT") {
    n <- ncol(res$gt)
    res$gt <- res$gt[, seq(1, n, 2)] + res$gt[, seq(2, n, 2)]
    res$gt[res$gt < 0] <- NA
  } else if (format == "GT") {
    res$gt[res$gt < 0] <- NA
  }
  return(res)
}
