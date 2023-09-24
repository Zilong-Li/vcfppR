#' @title
#' read VCF/BCF contents into R data structure
#'
#' @description
#' The swiss army knife to read VCF/BCF into R data types rapidly and easily.
#'
#' @details
#' vcftable("input.vcf.gz", "chr21:1-60000000", "id01, id02")
#' the similar interface in bcftools is =>
#' bcftools view -s "id01,id02" input.bcf.gz chr1:100000-20000
#'
#' @param vcffile path to the VCF/BCF file
#'
#' @param region region to subset like bcftools
#'
#' @param samples samples to subset like bcftools
#'
#' @param vartype restrict to specific type of variants, "snps","indels", "multisnps","multiallelics"
#'
#' @param format the FORMAT tag to extract, valid values are
#'               "GT","GP", "DP","DS","GL","PL","GQ","HQ","MQ","PQ"
#'
#' @param ploidy the ploidy of the organism if it is not diploidy
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
vcftable <- function(vcffile, region = "", samples = "-", vartype = "all", format = "GT", ploidy = 2) {
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
                DS = tableDS(vcffile, region, samples, snps, indels, multiallelics, multisnps),
                DP = tableDP(vcffile, region, samples, snps, indels, multiallelics, multisnps),
                PL = tablePL(vcffile, region, samples, snps, indels, multiallelics, multisnps),
                GL = tableGL(vcffile, region, samples, snps, indels, multiallelics, multisnps),
                HQ = tableHQ(vcffile, region, samples, snps, indels, multiallelics, multisnps),
                MQ = tableMQ(vcffile, region, samples, snps, indels, multiallelics, multisnps),
                PQ = tablePQ(vcffile, region, samples, snps, indels, multiallelics, multisnps),
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
