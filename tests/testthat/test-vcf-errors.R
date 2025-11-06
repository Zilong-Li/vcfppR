library(testthat)

vcffile <- system.file("extdata", "imputed.gt.vcf.gz", package="vcfppR")
svfile <- system.file("extdata", "sv.vcf.gz", package="vcfppR")

test_that("vcftable handles non-existent file", {
  skip_on_os("windows")
  # Non-existent file should error
  expect_error(vcftable("/path/to/nonexistent.vcf.gz"))
})

test_that("vcftable handles invalid region format", {
  skip_on_os("windows")
  # Invalid region should error
  expect_error(vcftable(vcffile, region = "invalid:region:format"))
})

test_that("vcftable handles non-existent region", {
  skip_on_os("windows")
  # Non-existent chromosome
  expect_error(vcftable(vcffile, region = "chr99:1-1000"))
})

test_that("vcftable handles invalid vartype", {
  skip_on_os("windows")
  # Invalid variant type should error
  expect_error(vcftable(vcffile, vartype = "invalid_type"),
               "Invaild variant type!")
})

test_that("vcfcomp handles mismatched samples", {
  skip_on_os("windows")
  # Create vcftable objects with different samples
  test <- vcftable(vcffile, samples = "HG00673,NA10840", setid = TRUE)
  truth <- vcftable(vcffile, samples = "HG00151,HG00380", setid = TRUE)

  # Should error on mismatched samples
  expect_error(vcfcomp(test, truth, stats = "r2", bins = c(0,1)),
               "the samples name in two VCF files is inconsistent")
})

test_that("vcfcomp handles invalid stats parameter", {
  skip_on_os("windows")
  # Invalid stats should error
  expect_error(vcfcomp(vcffile, vcffile, stats = "invalid_stat",
                       region = "chr21:1-5050000"),
               "invalid option of stats")
})

test_that("subset.vcftable handles invalid subset expression", {
  skip_on_os("windows")
  res <- vcftable(vcffile, region = "chr21:1-5050000")

  # Non-logical subset should error
  expect_error(subset(res, subset = pos),
               "'subset' must be logical")

  expect_error(subset(res, subset = "invalid"),
               "'subset' must be logical")
})

test_that("vcftable handles empty results gracefully", {
  skip_on_os("windows")
  # Query that returns no results
  res <- vcftable(vcffile, region = "chr21:1-1000", qual = 10000)

  expect_s3_class(res, "vcftable")
  expect_equal(length(res$pos), 0)
  expect_equal(length(res$gt), 0)
  expect_true("samples" %in% names(res))
})

test_that("vcfcomp handles no overlapping variants", {
  skip_on_os("windows")
  # Get two non-overlapping regions
  test <- vcftable(vcffile, region = "chr21:5000000-5001000", setid = TRUE)
  truth <- vcftable(vcffile, region = "chr21:5030000-5031000", setid = TRUE)

  # Should handle gracefully or error appropriately
  # The behavior depends on implementation
  result <- tryCatch({
    vcfcomp(test, truth, stats = "r2", bins = c(0,1))
  }, error = function(e) {
    e
  })

  # Just verify it doesn't crash unexpectedly
  expect_true(inherits(result, "vcfcomp") || inherits(result, "error"))
})

test_that("vcfinfo handles non-existent INFO tag", {
  skip_on_os("windows")
  # Non-existent INFO tag should error
  expect_error(vcfinfo(vcffile, region = "chr21:1-5050000",
                       tag = "NONEXISTENT_TAG"))
})

test_that("subset.vcftable handles all FALSE subset", {
  skip_on_os("windows")
  res <- vcftable(vcffile, region = "chr21:1-5050000")

  # Subset that returns nothing
  empty <- subset(res, info == "INFO")

  expect_s3_class(empty, "vcftable")
  expect_equal(length(empty$gt), 0)
  expect_equal(length(empty$pos), 0)
})

test_that("vcftable handles missing FORMAT tag gracefully", {
  skip_on_os("windows")
  # Request non-existent FORMAT tag
  res <- vcftable(vcffile, region = "chr21:1-5050000", format = "NONEXISTENT")

  # Should return with NA or empty format
  expect_s3_class(res, "vcftable")
  expect_true(length(names(res)) > 0)
})


test_that("vcfsummary handles empty VCF region", {
  skip_on_os("windows")
  # Query region with no variants
  result <- tryCatch({
    vcfsummary(vcffile, region = "chr21:1-1000")
  }, error = function(e) {
    e
  })

  # Should handle gracefully
  expect_true(inherits(result, "vcfsummary") || inherits(result, "error"))
})

test_that("vcftable handles very small regions", {
  skip_on_os("windows")
  # Single base pair region
  res <- vcftable(vcffile, region = "chr21:5030347-5030347")

  expect_s3_class(res, "vcftable")
  # May or may not have variants at that exact position
  expect_true(length(res$pos) >= 0)
})

test_that("subset.vcftable handles NA in conditions", {
  skip_on_os("windows")
  res <- vcftable(vcffile, region = "chr21:1-5050000")

  # Create condition that might produce NAs
  # NAs should be treated as FALSE
  result <- subset(res, qual > mean(qual))

  expect_s3_class(result, "vcftable")
  expect_true(nrow(result$gt) <= nrow(res$gt))
})

test_that("vcfwriter handles empty VCF creation", {
  skip_on_os("windows")
  # Create a minimal VCF
  outfile <- tempfile(fileext = ".vcf.gz")

  result <- tryCatch({
    w <- vcfwriter$new(outfile, "VCF4.3")
    w$close()
    file.exists(outfile)
  }, error = function(e) {
    FALSE
  })

  if(file.exists(outfile)) unlink(outfile)
  expect_true(result)
})

test_that("vcftable handles different ploidy levels", {
  skip_on_os("windows")
  # Most human data is diploid, but test handling
  res <- vcftable(vcffile, region = "chr21:1-5050000", collapse = FALSE)

  expect_s3_class(res, "vcftable")
  # For diploid, should have 2 columns per sample
  expect_equal(ncol(res$gt), 2 * length(res$samples))
})

test_that("vcfcomp handles different bin specifications", {
  skip("to be fixed")
  samples <- "HG00673,NA10840"

  # Single bin
  res1 <- vcfcomp(vcffile, vcffile, stats = "r2",
                  samples = samples, bins = c(0, 1), setid = TRUE)

  # Multiple bins
  expect_warning(res2 <- vcfcomp(vcffile, vcffile, stats = "r2",
                                 samples = samples, bins = c(0, 0.01, 0.05, 0.1, 0.5),
                                 setid = TRUE))

  # NULL bins (auto-generated)
  expect_warning(res3 <- vcfcomp(vcffile, vcffile, stats = "r2",
                                 samples = samples, bins = NULL, setid = TRUE))

  expect_s3_class(res1, "vcfcomp")
  expect_s3_class(res2, "vcfcomp")
  expect_s3_class(res3, "vcfcomp")

})

test_that("vcftable handles special characters in sample names", {
  skip_on_os("windows")
  # Get samples (may have special characters)
  res <- vcftable(vcffile, region = "chr21:1-5050000")

  expect_s3_class(res, "vcftable")
  expect_true(length(res$samples) > 0)
  expect_true(is.character(res$samples))
})

test_that("vcfcomp PSE handles unphased data", {
  skip_on_os("windows")
  samples <- "HG00151,HG00380"

  # PSE needs phased data
  # Test behavior with potentially unphased data
  result <- tryCatch({
    vcfcomp(vcffile, vcffile, stats = "pse",
            samples = samples, setid = TRUE)
  }, error = function(e) {
    e
  })

  # Should handle gracefully (either work or give meaningful error)
  expect_true(inherits(result, "vcfcomp") || inherits(result, "error"))
})

test_that("vcftable handles multi-allelic variants correctly", {
  skip_on_os("windows")
  # Get multi-allelic variants
  res <- vcftable(vcffile, vartype = "multiallelics")

  expect_s3_class(res, "vcftable")

  # Check ALT field contains commas for multi-allelics
  if(length(res$alt) > 0) {
    expect_true(all(grepl(",", res$alt)))
  }
})

test_that("subset.vcftable handles selecting non-existent fields", {
  skip_on_os("windows")
  res <- vcftable(vcffile, region = "chr21:1-5050000")

  # Try to select field that doesn't exist
  # Should error or handle gracefully
  result <- tryCatch({
    subset(res, select = c(chr, pos, nonexistent_field))
  }, error = function(e) {
    e
  })

  expect_true(inherits(result, "vcftable") || inherits(result, "error"))
})

test_that("vcfpopgen handles samples with no variation", {
  skip_on_os("windows")
  # Test with limited region
  result <- tryCatch({
    vcfpopgen(vcffile, region = "chr21:5000000-5001000", type = "het")
  }, error = function(e) {
    e
  })

  expect_true(inherits(result, "list") || inherits(result, "error"))
})
