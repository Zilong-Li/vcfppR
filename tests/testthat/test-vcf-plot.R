library(testthat)

vcffile <- system.file("extdata", "raw.gt.vcf.gz", package="vcfppR")
svfile <- system.file("extdata", "sv.vcf.gz", package="vcfppR")

test_that("vcfplot rejects invalid object types", {
  skip_on_os("windows")

  # Test with invalid object
  expect_error(vcfplot(list(a = 1, b = 2)),
               "the input object is not vcftable, vcfcomp or vcfsummary class")
  expect_error(vcfplot(data.frame(x = 1:10)),
               "the input object is not vcftable, vcfcomp or vcfsummary class")
})

test_that("vcfplot works with vcftable objects", {
  skip_on_os("windows")

  res <- vcftable(vcffile, "chr21:1-5050000", format = "GQ")

  # Should not error with default parameters
  expect_silent(vcfplot(res, which.sample = 1))

  # Should work with different samples
  expect_silent(vcfplot(res, which.sample = 2))
})

test_that("vcfplot with vcftable handles NULL which.sample", {
  skip_on_os("windows")

  res <- vcftable(vcffile, "chr21:1-5050000", format = "DP")

  # Should default to sample 1 and show message
  expect_message(vcfplot(res, which.sample = NULL),
                 "which.sample is NULL. will set which.sample as 1")
})

test_that("vcfplot with vcftable handles multiple samples in which.sample", {
  skip_on_os("windows")

  res <- vcftable(vcffile, "chr21:1-5050000", format = "GQ")

  # Should select first sample and show message
  expect_message(vcfplot(res, which.sample = c(1, 2)),
                 "which.sample has length > 1. select the first one")
})

test_that("vcfplot with vcftable accepts custom parameters", {
  skip_on_os("windows")

  res <- vcftable(vcffile, "chr21:1-5050000", format = "GQ")

  # Should work with custom plotting parameters
  expect_silent(vcfplot(res, which.sample = 1,
                        col = "red",
                        pch = 19,
                        main = "Custom Title",
                        xlab = "Position",
                        ylab = "Quality"))
})

test_that("vcfplot works with vcfsummary objects", {
  skip_on_os("windows")

  # Create summary
  res <- vcfsummary(vcffile, "chr21:1-5050000")

  # Should work without pop file (barplot)
  expect_silent(vcfplot(res))
})

test_that("vcfplot with vcfsummary accepts custom barplot parameters", {
  skip_on_os("windows")

  res <- vcfsummary(vcffile, "chr21:1-5050000")

  # Should work with custom parameters
  expect_silent(vcfplot(res,
                        col = "blue",
                        main = "Variant Summary",
                        xlab = "Type",
                        ylab = "Count"))
})

test_that("vcfplot works with vcfcomp r2 objects", {
  skip_on_os("windows")

  # Create a comparison
  comp <- vcfcomp(vcffile, vcffile,
                  stats = "r2", formats = c('GT', "GT" ),
                  region = "chr21:1-5050000")

  # Should work with default parameters
  expect_silent(vcfplot(comp))
})


test_that("vcfplot works with vcfcomp nrc objects", {
  skip_on_os("windows")

  comp <- vcfcomp(vcffile, vcffile,
                  stats = "nrc",
                  region = "chr21:1-5050000")

  expect_silent(vcfplot(comp))
})



test_that("vcfplot works with vcftable GT format", {
  skip_on_os("windows")

  res <- vcftable(vcffile, "chr21:1-5050000", format = "GT")

  # GT is collapsed by default, should work
  expect_silent(vcfplot(res, which.sample = 1))
})



test_that("vcfplot works with structural variants", {
  skip_on_os("windows")

  # Summary of SVs
  sv_summary <- vcfsummary(svfile, svtype = TRUE)

  # Should plot SV summary
  expect_silent(vcfplot(sv_summary))
})

test_that("vcfplot with vcftable works with different FORMAT tags", {
  skip_on_os("windows")

  # Test different FORMAT tags
  res_gq <- vcftable(vcffile, "chr21:1-5050000", format = "GQ")
  res_dp <- vcftable(vcffile, "chr21:1-5050000", format = "DP")
  res_pl <- vcftable(vcffile, "chr21:1-5050000", format = "PL", collapse = TRUE)

  expect_silent(vcfplot(res_gq, which.sample = 1))
  expect_silent(vcfplot(res_dp, which.sample = 1))
  expect_silent(vcfplot(res_pl, which.sample = 1))
})

test_that("vcfplot comp handles variant-wise statistics", {
  skip_on_os("windows")

  # by.variant = TRUE should produce different structure
  comp <- vcfcomp(vcffile, vcffile,
                  stats = "r2", formats = c("GT", "GT" ),
                  by.variant = TRUE,
                  region = "chr21:1-5050000")

  # Should work without specifying which.sample
  expect_silent(vcfplot(comp))
})


## vcffiles <- rep("/home//zilong/Downloads/longcallD.vcf.gz", 5)
## region <- "chr1:570000-592000"
## plot_variants_per_haplotype(vcffiles, region)



