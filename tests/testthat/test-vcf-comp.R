library(testthat)

rawvcf <- system.file("extdata", "raw.gt.vcf.gz", package="vcfppR")
imputedvcf <- system.file("extdata", "imputed.gt.vcf.gz", package="vcfppR")

test_that("can work for correlation r2 between DS and GT", {
  skip_on_os(c("windows"), arch = NULL)
  samples <- "HG00673,NA10840"
  res <- vcfcomp(imputedvcf, imputedvcf, bins = c(0,1), by.sample = TRUE, samples = samples)
  expect_equal(as.numeric(res$r2[1,]), c(15, 6, 1), tolerance=1e-6)
  expect_identical(unlist(res$f1[[3]]), rep(1,2))
  expect_identical(unlist(res$nrc[[3]]), rep(1,2))
})

test_that("can work for F1 score", {
  skip_on_os(c("windows"), arch = NULL)
  samples <- "HG00673,NA10840"
  res <- vcfcomp(rawvcf, imputedvcf, stats = "f1", samples = samples, by.sample = TRUE, bins = c(0,1))
  expect_identical(paste0(res$samples, collapse = ","), samples)  
  expect_identical(unlist(res[[2]][[3]]), rep(1,2))
})

test_that("can work for NRC rate by sample", {
  skip_on_os(c("windows"), arch = NULL)
  samples <- "HG00673,NA10840"
  res <- vcfcomp(rawvcf, imputedvcf, stats = "nrc", samples = samples, by.sample = TRUE, bins = c(0,1))
  expect_identical(paste0(res$samples, collapse = ","), samples)  
  expect_identical(unlist(res[[2]][[3]]), rep(1,2))
})

test_that("can work for PSE", {
  skip_on_os(c("windows"), arch = NULL)
  samples <- "HG00151,HG00380"
  expect_warning(res <- vcfcomp(imputedvcf, imputedvcf, stats = "pse", samples = samples))
  expect_true(is.nan(res[[2]][[1]]$pse))
  expect_true(is.na(res[[2]][[2]]))
  expect_identical(paste0(res$samples, collapse = ","), samples)
})


