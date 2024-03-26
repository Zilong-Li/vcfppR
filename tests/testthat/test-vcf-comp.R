library(testthat)

rawvcf <- system.file("extdata", "raw.gt.vcf.gz", package="vcfppR")
imputedvcf <- system.file("extdata", "imputed.gt.vcf.gz", package="vcfppR")

test_that("can work for correlation r2 between DS and GT", {
  skip_on_os(c("windows"), arch = NULL)
  res <- vcfcomp(imputedvcf, imputedvcf)
  expect_identical(as.numeric(res$r2[5,]), c(2, 2, 1))
  expect_equal(as.numeric(res$r2[8,]), c(5, 49, 0.9838495), tolerance=1e-6)
})

test_that("can work for F1 score", {
  skip_on_os(c("windows"), arch = NULL)
  samples <- "HG00673,NA10840"
  res <- vcfcomp(rawvcf, imputedvcf, stats = "f1", samples = samples, bysample = T, bins = c(0,1))
  expect_identical(paste0(res$samples, collapse = ","), samples)  
  expect_identical(unlist(res$f1[[3]]), rep(1,2))
})

test_that("can work for NRC rate by sample", {
  skip_on_os(c("windows"), arch = NULL)
  samples <- "HG00673,NA10840"
  res <- vcfcomp(rawvcf, imputedvcf, stats = "nrc", samples = samples, bysample = T, bins = c(0,1))
  expect_identical(paste0(res$samples, collapse = ","), samples)  
  expect_identical(unlist(res$nrc[[3]]), rep(1,2))
})

