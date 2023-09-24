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
#' @examples
#' library('vcfppR')
#' vcffile <- "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr21.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
#' res <- popgen.heterozygosity(vcffile, "chr21:1-10000000")
#' str(res)
#' ped <- read.table("https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt", h=T)
#' ped <- ped[order(ped$Superpopulation),]
#' out <- sapply(unique(ped$Superpopulation), function(pop) {
#'   id <- subset(ped, Superpopulation == pop)[,"SampleID"]
#'   ord <- match(id, res$samples)
#'   res$hets[ord]
#' })
#' boxplot(out)
#' @export
popgen.heterozygosity <- function(vcffile, region = "", samples = "-", pass = FALSE, qual = 0) {
    return(heterozygosity(vcffile, region, samples, pass, qual))
}
