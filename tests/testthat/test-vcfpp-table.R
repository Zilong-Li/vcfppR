library(vcfppR)


vcffile <- system.file("extdata", "test-GL.vcf.gz", package="vcfppR")

(res <- vcftable(vcffile,"chr20", info = FALSE, format = "na"))
str(res)

(res <- vcftable(vcffile,"chr20"))

vcffile <- system.file("extdata", "test-PL.vcf.gz", package="vcfppR")
(res <- vcftable(vcffile,"chr20", format = "PL"))

vcffile <- "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr21.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
res <- vcftable(vcffile, "chr21:1-5100000", vartype = "snps")

q()

vcffile <- "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr21.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
res <- popgen.heterozygosity(vcffile, "chr21:1-10000000")
str(res)


ped <- read.table("https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt", h=T)
ped <- ped[order(ped$Superpopulation),]

out <- sapply(unique(ped$Superpopulation), function(pop) {
  id <- subset(ped, Superpopulation == pop)[,"SampleID"]
  ord <- match(id, res$samples)
  res$hets[ord]
})
boxplot(out, main = "Heterozygous variants")


vcffile <- "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr21.recalibrated_variants.vcf.gz"
ped <- read.table("https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt", h=T)
ped <- ped[order(ped$Superpopulation),]

region <- "chr21:10000000-10500000"
res <- vcfsummary(vcffile, region)
str(res)
out <- sapply(unique(ped$Superpopulation), function(pop) {
  id <- subset(ped, Superpopulation == pop)[,"SampleID"]
  ord <- match(id, res$samples)
  res$SNP[ord] + res$INDEL[ord]
})

## devtools::build_readme()

boxplot(out, main = paste0("Average number of SNP & INDEL variants\nin region ", region))

dev.off()

svfile <- "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20210124.SV_Illumina_Integration/1KGP_3202.gatksv_svtools_novelins.freeze_V3.wAF.vcf.gz"
sv <- vcfsummary(svfile, svtype = TRUE)
str(sv)

allsvs <- sv$summary[-1]
bar <- barplot(allsvs, ylim = c(0, 1.1*max(allsvs)), main = "Variant Count (all SVs)")
text(bar, allsvs+4500, paste("n: ", allsvs, sep=""))

save.image("~/tmp/vcfppR.RData")

load("~/tmp/vcfppR.RData")
ls()

samples <- sv$samples
supers <- unique(ped$Superpopulation)
o <- lapply(supers, function(pop) {
  id <- subset(ped, Superpopulation == pop)[,"SampleID"]
  ord <- match(id, samples)
  sapply(sv[-c(1,2)], "[", ord)
})
names(o) <- supers

## s <- as.data.frame(lapply(o, colSums))
s <- as.data.frame(lapply(o, colMeans))
barplot(colSums(s))

## s[, c("AFR", "EUR", "SAS", "EAS", "AMR")]

s <- as.data.frame(formatC(as.matrix(s), digits = 2, format = "fg"))

barplot(as.matrix(s))

dev.off()

