library(testthat)

rawvcf <- system.file("extdata", "raw.gt.vcf.gz", package="vcfppR")
imputedvcf <- system.file("extdata", "imputed.gt.vcf.gz", package="vcfppR")

test_that("can work for F1 score", {
  skip_on_os(c("windows"), arch = NULL)
  samples <- "HG00673,NA10840"
  res <- vcfcomp(rawvcf, imputedvcf, stats = "f1", samples = samples)
  expect_identical(paste0(res$samples, collapse = ","), samples)  
  expect_identical(res$f1, rep(1,2))
})

test_that("can work for correlation r2 between DS and GT", {
  skip_on_os(c("windows"), arch = NULL)
  res <- vcfcomp(imputedvcf, imputedvcf)
  res$r2[8,]
  rownames(res$r2)
  expect_identical(as.numeric(res$r2[5,]), c(2, 2, 1, 1))
  expect_equal(as.numeric(res$r2[8,]), c(5, 49, 0.9838495, 0.9838479), tolerance=1e-6)
})
