library(testthat)

vcffile <- system.file("extdata", "raw.gt.vcf.gz", package="vcfppR")

test_that("vcfreader works", {
  br <- vcfreader$new(vcffile, "", "")
  str(br)
  br$getvariant()
  gt <- br$genotypes()
  expect_identical(length(gt), 6404L)
})
