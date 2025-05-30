library(testthat)

rawvcf <- system.file("extdata", "raw.gt.vcf.gz", package="vcfppR")
imputedvcf <- system.file("extdata", "imputed.gt.vcf.gz", package="vcfppR")

test_that("can work for correlation r2 between DS and GT", {
  skip_on_os(c("windows"), arch = NULL)
  samples <- "HG00673,NA10840"
  res <- vcfcomp(imputedvcf, imputedvcf, bins = c(0,1),
                 by.sample = TRUE, samples = samples, stats = "all",
                 setid = TRUE)
  expect_equal(as.numeric(unlist(res$r2[1,])), c(15, 6, 1, 1), tolerance=1e-6)
  expect_identical(unlist(res$f1[[3]]), rep(1,2))
  expect_identical(unlist(res$nrc[[3]]), rep(1,2))
})

test_that("can work for r2 with af", {
  skip_on_os(c("windows"), arch = NULL)
  samples <- "HG00673,NA10840"
  d1 <- vcftable(imputedvcf,setid = T, info=F)
  af <- runif(15)
  names(af) <- d1$id
  affile <- tempfile()
  saveRDS(af, affile)
  res <- vcfcomp(imputedvcf, rawvcf, stats = "r2", bins = c(0,1),
                 af = affile, samples = samples, setid = TRUE)
  expect_identical(res[[2]][,3], 1)
})

test_that("can work for F1 score", {
  skip_on_os(c("windows"), arch = NULL)
  samples <- "HG00673,NA10840"
  res <- vcfcomp(rawvcf, imputedvcf, stats = "f1", samples = samples,
                 by.sample = TRUE, bins = c(0,1), setid = TRUE)
  expect_identical(paste0(res$samples, collapse = ","), samples)  
  expect_identical(unlist(res[[2]][[3]]), rep(1,2))
})

test_that("can work for NRC rate by sample", {
  skip_on_os(c("windows"), arch = NULL)
  samples <- "HG00673,NA10840"
  res <- vcfcomp(rawvcf, imputedvcf, stats = "nrc", samples = samples,
                 by.sample = TRUE, bins = c(0,1), setid = TRUE)
  expect_identical(paste0(res$samples, collapse = ","), samples)  
  expect_identical(unlist(res[[2]][[3]]), rep(1,2))
})

test_that("can work for PSE", {
  skip_on_os(c("windows"), arch = NULL)
  samples <- "HG00151,HG00380"
  (res <- vcfcomp(imputedvcf, imputedvcf, stats = "pse", samples = samples, setid = TRUE))
  expect_true(is.nan(res[[2]][[1]]$pse))
  expect_true((res[[2]][[2]]$pse==0))
  expect_true((res[[2]][[2]]$disc==0))
  expect_true(is.null(res[[2]][[2]]$pos))
  expect_identical(paste0(res$samples, collapse = ","), samples)
})



test_that("can work for vcftable objects", {
  skip_on_os(c("windows"), arch = NULL)
  samples <- "HG00673,NA10840"
  test <- vcftable(rawvcf, samples = samples, setid = T)
  truth <- vcftable(imputedvcf, samples = samples, setid = T)
  res <- vcfcomp(test, truth, stats = "nrc", samples = samples,
                 by.sample = TRUE, bins = c(0,1), setid = TRUE)
  expect_identical(paste0(res$samples, collapse = ","), samples)  
  expect_identical(unlist(res[[2]][[3]]), rep(1,2))
})
