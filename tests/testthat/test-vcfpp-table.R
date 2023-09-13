library(vcfppR)

vcffile <- system.file("extdata", "test-GL.vcf.gz", package="vcfppR")
(gt <- tableGT(vcffile,"chr20"))
(gl <- tableGL(vcffile,"chr20"))

vcffile <- system.file("extdata", "test-PL.vcf.gz", package="vcfppR")
(pl <- tablePL(vcffile,"chr20"))

(hets <- heterozygosity(vcffile))



## vcffile <- "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr21.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
vcffile <- "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr21.recalibrated_variants.vcf.gz"
res <- summaryVariants(vcffile, "chr21:1-10000000")

do.call(cbind, res)
