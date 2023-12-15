library(vcfppR)

test_that("summarize SVs", {
  svfile <- system.file("extdata", "sv.vcf.gz", package="vcfppR")
  res <- vcfsummary(svfile, svtype = TRUE)
  expect_equal(as.integer(res$summary["SV"]), 144L)
  expect_equal(as.integer(res$summary["DUP"]), 36L)
})


test_that("summarize small variants", {
  rawvcf <- system.file("extdata", "raw.gt.vcf.gz", package="vcfppR")
  res <- vcfsummary(rawvcf)
  expect_equal(as.integer(res$summary["ALL"]), 72L)
  expect_equal(as.integer(res$summary["SNP"]), 63L)
  expect_equal(as.integer(res$summary["INDEL"]), 5L)
  expect_equal(as.integer(res$summary["MultiAllelics"]), 4L)
  expect_equal(as.integer(res$summary["MultiAllelicSNP"]), 3L)
})



