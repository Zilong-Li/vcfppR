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
#'
#' @param format the FORMAT tag to extract, valid values are
#'               "GT","GP", "DP","DS","GL","PL","GQ","HQ","MQ","PQ"
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
#' vcffile <- system.file("extdata", "test-GL.vcf.gz", package="vcfppR")
#' (res <- vcftable(vcffile,"chr20"))
#' vcffile <- system.file("extdata", "test-PL.vcf.gz", package="vcfppR")
#' (res <- vcftable(vcffile,"chr20", format = "PL"))
#' vcffile <- "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr21.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
#' res <- vcftable(vcffile, "chr21:1-5100000", vartype = "snps")
#' str(res)
#' @export
vcftable <- function(vcffile, region = "", samples = "-", vartype = "all", format = "GT") {
  snps <- FALSE
  indels <- FALSE
  multiallelics <- FALSE
  multisnps <- FALSE
  if(vartype == "snps") snps <- TRUE
  else if(vartype == "indels") indels <- TRUE
  else if(vartype == "multisnps") multisnps <- TRUE
  else if(vartype == "multiallelics") multiallelics <- TRUE
  else if(vartype != "all") stop("Invaild variant type!")
  res <- switch(format,
                GT = tableGT(vcffile, region, samples, snps, indels, multiallelics, multisnps),
                GP = tableGP(vcffile, region, samples, snps, indels, multiallelics, multisnps),
                GQ = tableGQ(vcffile, region, samples, snps, indels, multiallelics, multisnps),
                AD = tableAD(vcffile, region, samples, snps, indels, multiallelics, multisnps),
                DS = tableDS(vcffile, region, samples, snps, indels, multiallelics, multisnps),
                DP = tableDP(vcffile, region, samples, snps, indels, multiallelics, multisnps),
                PL = tablePL(vcffile, region, samples, snps, indels, multiallelics, multisnps),
                GL = tableGL(vcffile, region, samples, snps, indels, multiallelics, multisnps),
                HQ = tableHQ(vcffile, region, samples, snps, indels, multiallelics, multisnps),
                MQ = tableMQ(vcffile, region, samples, snps, indels, multiallelics, multisnps),
                PQ = tablePQ(vcffile, region, samples, snps, indels, multiallelics, multisnps),
                stop("Invaild tag in FORAMT column"))
  res[[10]] <- do.call("rbind", res[[10]])
  if(format == "GT") {
    n <- ncol(res$gt)
    ploidy <- n / length(res$samples)
    if(ploidy == 2) {
      res$gt <- res$gt[, seq(1, n, 2)] + res$gt[, seq(2, n, 2)]
      res$gt[res$gt < 0] <- NA
    } else {
      res$gt[res$gt < 0] <- NA
    }
  }
  return(res)
}
