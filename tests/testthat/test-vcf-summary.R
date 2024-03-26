library(vcfppR)

test_that("summarize SVs", {
  svfile <- system.file("extdata", "sv.vcf.gz", package="vcfppR")
  res <- vcfsummary(svfile, svtype = TRUE)
  expect_equal(as.integer(res$summary["SV"]), 9L)
  expect_equal(as.integer(res$summary["DUP"]), 4L)
})


test_that("summarize small variants", {
  rawvcf <- system.file("extdata", "raw.gt.vcf.gz", package="vcfppR")
  res <- vcfsummary(rawvcf)
  expect_equal(as.integer(res$summary["ALL"]), 15L)
  expect_equal(as.integer(res$summary["SNP"]), 12L)
  expect_equal(as.integer(res$summary["INDEL"]), 1L)
  expect_equal(as.integer(res$summary["MNP"]), 0L)
  expect_equal(as.integer(res$summary["MultiAllelics"]), 2L)
  expect_equal(as.integer(res$summary["MultiAllelicSNP"]), 2L)
})



