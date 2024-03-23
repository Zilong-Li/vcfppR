library(testthat)

test_that("vcfwriter: writing variant works", {
  ## skip_on_os(c("windows"), arch = NULL)
  outvcf <- file.path(paste0(tempfile(), ".vcf.gz"))
  bw <- vcfwriter$new(outvcf, "VCF4.3")
  bw$addLine('##ALT=<ID=NON_REF,Description="Represents any possible alternative allele at this location">')
  bw$addContig("chr20")
  bw$addFILTER("PASS", "All filters passed")
  bw$addINFO("AF", "A", "Float", "Estimated allele frequency in the range (0,1)");
  s1 <- "chr20\t2006060\t.\tG\tC\t100\tPASS\t."
  bw$writeline(s1)
  bw$close()
  ## tests
  br <- vcfreader$new(outvcf)
  print(br$header())
  br$variant()
  s2 <- gsub("\n", "", br$string())
  expect_identical(br$chr(), "chr20")
  expect_identical(s1, s2)
  s2 <- br$line()
  expect_identical(s1, s2)
})


