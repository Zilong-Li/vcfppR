library(testthat)

vcffile <- system.file("extdata", "raw.gt.vcf.gz", package="vcfppR")
svfile <- system.file("extdata", "sv.vcf.gz", package="vcfppR")

test_that("extract FORMAT item not exists", {
  res <- vcftable(vcffile, vartype = "snps", format = "AA")
  expect_identical(res$na, as.integer(NA))
})

test_that("extract GT for all SNPs", {
  res <- vcftable(vcffile, vartype = "snps")
  ## if all are SNPs and ALT!="."
  expect_identical(sum(res$alt!=""), length(res$alt))
})

test_that("extract SNPs by ID and reset ID afterwards", {
  res <- vcftable(vcffile, vartype = "snps", id = c("chr21:5030516:G:A"), setid = TRUE)
  expect_identical(res$id, "chr21_5030516_G_A")
})

test_that("extract GT for a indel with FILTER not displaying PASS", {
  res <- vcftable(vcffile, id = c("chr21:5030240:AC:A"), vartype = "indels", pass = TRUE)
  expect_equal(length(res$gt), 0)
})

test_that("extract GT for variant with ID=chr21:5030516:G:A", {
  res <- vcftable(vcffile, id = c("chr21:5030516:G:A"), vartype = "snps")
  expect_equal(res$qual, 222.25)
  expect_equal(res$ref, "G")
  expect_equal(res$alt, "A")
})

test_that("extract GT for all multisnps", {
  res <- vcftable(vcffile, vartype = "multisnps")
  expect_identical(res$ref, c("C", "C"))
  expect_identical(res$alt, c("G,T", "T,G"))
})

test_that("extract GT for all multiallelics", {
  res <- vcftable(vcffile, vartype = "multiallelics")
  expect_identical(res$ref, c("C", "C"))
  expect_identical(res$alt, c("G,T", "T,G"))
})

test_that("extract AD (Integer, number=.) for variant with ID=chr21:5030247:G:A", {
  res <- vcftable(vcffile, id = c("chr21:5030347:G:A"), pass = TRUE, format = "AD")
  expect_equal(nrow(res$AD), 1)
  expect_equal(ncol(res$AD), length(res$samples) * 2L)
  expect_equal(res$AD[1,1:2], c(5, 0))
})

test_that("extract PL (Integer, number=G) for variant with ID=chr21:5030247:G:A", {
  res <- vcftable(vcffile, id = c("chr21:5030347:G:A"), pass = TRUE, format = "PL")
  expect_equal(nrow(res$PL), 1)
  expect_equal(ncol(res$PL), length(res$samples) * 3L)
  expect_equal(res$PL[1,1:3], c(0, 15, 174))
})

test_that("extract AB (Float, number=1)", {
  res <- vcftable(vcffile, id = c("chr21:5030347:G:A", "chr21:5030356:G:C"), format = "AB")
  expect_equal(nrow(res$AB), 2)
  expect_equal(ncol(res$AB), length(res$samples))
  expect_equal(all(is.na(res$AB[1,])), FALSE)
  expect_equal(all(is.na(res$AB[2,])), FALSE)
  expect_equal(res$AB[1,2893], 0.71, tolerance = 1e-6)
})

test_that("extract EV (String, number=1) from sv file", {
  ## res <- vcftable(svfile, format = "EV") ## there is no EV for ID=HGSV_240945
  res <- vcftable(svfile, vartype = "sv", format = "EV", id = c("HGSV_240934"))
  expect_equal(res$EV[1,], rep("RD", length(res$samples)))
})
