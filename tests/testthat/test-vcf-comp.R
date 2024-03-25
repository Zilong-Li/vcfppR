library(testthat)

rawvcf <- system.file("extdata", "raw.gt.vcf.gz", package="vcfppR")
imputedvcf <- system.file("extdata", "imputed.gt.vcf.gz", package="vcfppR")

test_that("can work for F1 score", {
  skip_on_os(c("windows"), arch = NULL)
  samples <- "HG00133,HG00262"
  res <- vcfcomp(rawvcf, imputedvcf, stats = "f1", format=c("GT", 'GT'), samples = samples)
  expect_identical(paste0(res$samples, collapse = ","), samples)  
  expect_identical(res$f1, rep(1,2))
})

test_that("can work for correlation r2 between DS and GT", {
  skip_on_os(c("windows"), arch = NULL)
  res <- vcfcomp(imputedvcf, imputedvcf)
  expect_identical(as.numeric(res$r2[5,]), c(16, 16, 1, 1))
})
