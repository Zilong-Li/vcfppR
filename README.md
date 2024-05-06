
<!-- README.md is generated from README.Rmd. Please edit that file -->

# vcfppR: rapid manipulation of the VCF/BCF file

<!-- badges: start -->

<a href="https://github.com/Zilong-Li/random/blob/main/vcfppR.png"><img src="https://raw.githubusercontent.com/Zilong-Li/random/main/vcfppR.png" width="120" align="right" /></a>
![R-CMD-check](https://github.com/Zilong-Li/vcfppR/actions/workflows/check-release.yaml/badge.svg)
[![CRAN
status](https://www.r-pkg.org/badges/version/vcfppR)](https://CRAN.R-project.org/package=vcfppR)
![<https://github.com/Zilong-Li/vcfppR/releases/latest>](https://img.shields.io/github/v/release/Zilong-Li/vcfppR.svg)
[![codecov](https://codecov.io/github/Zilong-Li/vcfppR/graph/badge.svg?token=QE1UFVRH98)](https://app.codecov.io/github/Zilong-Li/vcfppR)
[![Downloads](https://cranlogs.r-pkg.org/badges/vcfppR?color=blue)](https://CRAN.R-project.org/package=vcfppR)
<!-- badges: end -->

The vcfppR package implements various powerful functions for fast
genomics analyses with VCF/BCF files using the C++ API of
[vcfpp.h](https://github.com/Zilong-Li/vcfpp).

  - Load/save content of VCF/BCF into R objects with highly customizable
    options
  - Visualize and chracterize variants
  - Compare two VCF/BCF files and report various statistics
  - Streaming read/write VCF/BCF files with fine control of everything
  - [Paper](https://doi.org/10.1093/bioinformatics/btae049) shows vcfppR
    is 20x faster than `vcfR`. Also, much faster than `cyvcf2`

## Installation

``` r
## install.package("vcfppR") ## from CRAN
remotes::install_github("Zilong-Li/vcfppR") ## from latest github
```

## Cite the work

If you find it useful, please cite the
[paper](https://doi.org/10.1093/bioinformatics/btae049)

``` r
library(vcfppR)
citation("vcfppR")
```

## `vcftable`: read VCF as tabular data

`vcftable` gives you fine control over what you want to extract from
VCF/BCF files.

**Read only SNP variants**

``` r
library(vcfppR)
vcffile <- "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr21.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
res <- vcftable(vcffile, "chr21:1-5100000", vartype = "snps")
str(res)
```

**Read only SNP variants with PL format and drop the INFO column in the
VCF/BCF**

``` r
vcffile <- "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr21.recalibrated_variants.vcf.gz"
res <- vcftable(vcffile, "chr21:1-5100000", vartype = "snps", format = "PL", info = FALSE)
str(res)
```

**Read only INDEL variants with DP format in the VCF/BCF**

``` r
vcffile <- "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr21.recalibrated_variants.vcf.gz"
res <- vcftable(vcffile, "chr21:1-5100000", vartype = "indels", format = "DP")
str(res)
```

## `vcfcomp`: compare two VCF files and report concordance

Want to investigate the concordance between two VCF files? `vcfcomp` is
the utility function you need\!

**Genotype correlation**

``` r
vcffile <- "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr21.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
res <- vcfcomp(test = vcffile, truth = vcffile, region = "chr21:1-5100000", stats = "r2", formats = c('GT','GT'))
str(res)
```

**Genotype F1 score**

``` r
vcffile <- "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr21.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
res <- vcfcomp(test = vcffile, truth = vcffile, region = "chr21:1-5100000", stats = "f1")
str(res)
```

**Genotype Non-Reference Concordance**

``` r
vcffile <- "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr21.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
res <- vcfcomp(test = vcffile, truth = vcffile, region = "chr21:1-5100000", stats = "nrc")
str(res)
```

## `vcfsummary`: variants characterization

Want to summarize variants discovered by genotype caller e.g. GATK?
`vcfsummary` is the utility function you need\!

**Small variants**

``` r
vcffile <- "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr21.recalibrated_variants.vcf.gz"
region <- "chr21:10000000-10010000"
res <- vcfsummary(vcffile, region)
str(res)
# get labels and do plottiing
ped <- read.table("https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt", h=T)
ped <- ped[order(ped$Superpopulation),]
out <- sapply(unique(ped$Superpopulation), function(pop) {
  id <- subset(ped, Superpopulation == pop)[,"SampleID"]
  ord <- match(id, res$samples)
  res$SNP[ord] + res$INDEL[ord]
})

boxplot(out, main = paste0("Average number of SNP & INDEL variants\nin region ", region))
```

**Complex structure variants**

``` r
svfile <- "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20210124.SV_Illumina_Integration/1KGP_3202.gatksv_svtools_novelins.freeze_V3.wAF.vcf.gz"
sv <- vcfsummary(svfile, svtype = TRUE, region = "chr20")
str(sv)
allsvs <- sv$summary[-1]
bar <- barplot(allsvs, ylim = c(0, 1.1*max(allsvs)),
               main = "Variant Counts on chr20 (all SVs)")
```

## API

There are two classes i.e. ***vcfreader*** and ***vcfwriter*** offering
the full R-bindings of vcfpp.h. Check out the examples in the
[tests](tests/testthat) folder or refer to the manual.

``` r
?vcfppR::vcfreader
```
