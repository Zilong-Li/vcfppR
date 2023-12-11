library(vcfppR)
library(testthat)

vcffile <- system.file("extdata", "raw.gt.vcf.gz", package="vcfppR")

test_that("extract GT for all SNPs", {
  res <- vcftable(vcffile, vartype = "snps", format = "AA")
  ## if all are SNPs and ALT!="."
  expect_identical(sum(res$alt!=""), length(res$alt))
})

test_that("extract GT for variant with ID=chr21:5030516:G:A", {
  res <- vcftable(vcffile, id = c("chr21:5030516:G:A"))
  expect_equal(res$qual, 222.25)
  expect_equal(res$ref, "G")
  expect_equal(res$alt, "A")
})

