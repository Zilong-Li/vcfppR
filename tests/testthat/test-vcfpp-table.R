library(vcfppR)

vcffile <- system.file("extdata", "test-GL.vcf.gz", package="vcfppR")

(res <- vcftable(vcffile,"chr20"))
str(res)


vcffile <- system.file("extdata", "test-PL.vcf.gz", package="vcfppR")
(res <- vcftable(vcffile,"chr20", format = "PL"))

(hets <- heterozygosity(vcffile))



## vcffile <- "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr21.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
vcffile <- "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr21.recalibrated_variants.vcf.gz"

region <- "chr21:1-6030300"
region <- "chr21:10000000-11030300"
region <- "chr21:1-5130300"
res <- summaryVariants(vcffile, region)

str(res)
res$summary

svfile <- "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20210124.SV_Illumina_Integration/1KGP_3202.gatksv_svtools_novelins.freeze_V3.wAF.vcf.gz"
res <- summarySVs(svfile)
str(res)
res$summary

summary(res$CTX)

