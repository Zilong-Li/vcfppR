library(testthat)

rawvcf <- system.file("extdata", "raw.gt.vcf.gz", package="vcfppR")
imputedvcf <- system.file("extdata", "imputed.gt.vcf.gz", package="vcfppR")

test_that("can work for F1 score", {
  samples <- "HG00133,HG00143,HG00262"
  res <- vcfcomp(imputedvcf, imputedvcf, stats = "f1", format=c("GT", 'GT'), samples = samples)
  expect_identical(paste0(res$samples, collapse = ","), samples)  
  expect_identical(res$f1, rep(1,3))
})

test_that("can work for correlation r2 between DS and GT", {
  res <- vcfcomp(imputedvcf, imputedvcf)
  r2 <- res$r2
  expect_identical(as.numeric(res$r2[5,]), c(16, 16, 1, 1))
})
