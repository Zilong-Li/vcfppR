library(testthat)

bcffile <- system.file("extdata", "test.bcf", package="vcfppR")
vcffile <- system.file("extdata", "raw.gt.vcf.gz", package="vcfppR")
svfile <- system.file("extdata", "sv.vcf.gz", package="vcfppR")

test_that("extract INFO tag with Integer type", {
  res <- vcfinfo(vcffile, "AC", region = "chr21:1-5050000", vartype = 'snps')
  ## expect_identical(res$na, as.integer(NA))
  expect_equal(length(res$AC), 12L)
})

test_that("extract INFO tag with Float type", {
  res <- vcfinfo(vcffile, "AF", region = "chr21:1-5050000", vartype = 'snps')
  expect_equal(length(res$AF), 12L)
})

test_that("extract INFO tag with String type", {
  res <- vcfinfo(svfile, "SVTYPE", region = "chr21:1-5050000")
  expect_equal(length(res$SVTYPE), 2L)
})

test_that("extract INFO tag that does not exist", {
  expect_error(vcfinfo(svfile, "AA", region = "chr21:1-5050000"))
})

