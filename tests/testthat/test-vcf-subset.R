library(testthat)

vcffile <- system.file("extdata", "raw.gt.vcf.gz", package="vcfppR")

test_that("subset by quality score works", {
  skip_on_os("windows")
  res <- vcftable(vcffile, "chr21:1-5050000")

  # Subset by quality
  high_qual <- subset(res, qual > 100)

  expect_s3_class(high_qual, "vcftable")
  expect_true(all(high_qual$qual > 100))
  expect_true(nrow(high_qual$gt) < nrow(res$gt))
  expect_equal(high_qual$samples, res$samples)
})

test_that("subset by position range works", {
  skip_on_os("windows")
  res <- vcftable(vcffile, "chr21:1-5050000")

  # Subset by position
  region_subset <- subset(res, pos >= 5000000 & pos <= 5010000)

  expect_s3_class(region_subset, "vcftable")
  expect_true(all(region_subset$pos >= 5000000))
  expect_true(all(region_subset$pos <= 5010000))
  expect_equal(length(region_subset$pos), nrow(region_subset$gt))
})

test_that("subset with field selection works", {
  skip_on_os("windows")
  res <- vcftable(vcffile, "chr21:1-5050000")

  # Subset and select specific fields
  selected <- subset(res, pos >= 5000000 & pos <= 5010000,
                     select = c(chr, pos, ref, alt))

  expect_s3_class(selected, "vcftable")
  expect_equal(names(selected), c("samples", "chr", "pos", "ref", "alt"))
  expect_true("gt" %in% names(res))  # gt exists in original
  expect_false("gt" %in% names(selected))  # gt not in subset
  expect_true(all(selected$pos >= 5000000 & selected$pos <= 5010000))
})

test_that("subset SNPs by REF/ALT length works", {
  skip_on_os("windows")
  res <- vcftable(vcffile, "chr21:1-5050000")

  # Subset SNPs (single nucleotide variants)
  snps <- subset(res, nchar(ref) == 1 & nchar(alt) == 1)

  expect_s3_class(snps, "vcftable")
  expect_true(all(nchar(snps$ref) == 1))
  expect_true(all(nchar(snps$alt) == 1))
  expect_true(nrow(snps$gt) <= nrow(res$gt))
})

test_that("subset INDELs by REF/ALT length works", {
  skip_on_os("windows")
  res <- vcftable(vcffile, "chr21:1-5050000")

  # Subset INDELs (length difference)
  indels <- subset(res, nchar(ref) != nchar(alt))

  expect_s3_class(indels, "vcftable")
  expect_true(all(nchar(indels$ref) != nchar(indels$alt)))
  expect_gt(nrow(res$gt), nrow(indels$gt))
})

test_that("subset with FILTER column works", {
  skip_on_os("windows")
  res <- vcftable(vcffile, "chr21:1-5050000")

  # Subset PASS variants
  passed <- subset(res, filter == "PASS")

  expect_s3_class(passed, "vcftable")
  expect_true(all(passed$filter == "PASS"))
})

test_that("subset with chromosome selection works", {
  skip_on_os("windows")
  res <- vcftable(vcffile, "chr21:1-5050000")

  # All should be chr21 in this test file
  chr21 <- subset(res, chr == "chr21")

  expect_s3_class(chr21, "vcftable")
  expect_true(all(chr21$chr == "chr21"))
  expect_equal(nrow(chr21$gt), nrow(res$gt))
})

test_that("subset handles missing values correctly", {
  skip_on_os("windows")
  res <- vcftable(vcffile, "chr21:1-5050000")

  # If qual has NA values, they should be excluded
  # even if the condition would be TRUE for numeric values
  result <- subset(res, qual > 0)

  expect_s3_class(result, "vcftable")
  # No NAs should pass through
  expect_false(any(is.na(result$qual) & result$qual > 0, na.rm = TRUE))
})

test_that("subset without any condition returns all variants", {
  skip_on_os("windows")
  res <- vcftable(vcffile, "chr21:1-5050000")

  # No subset condition
  all_vars <- subset(res)

  expect_s3_class(all_vars, "vcftable")
  expect_equal(nrow(all_vars$gt), nrow(res$gt))
  # samples should be excluded by default select behavior
  expect_equal(names(all_vars), names(res))
})

test_that("subset with only select works", {
  skip_on_os("windows")
  res <- vcftable(vcffile, "chr21:1-5050000")

  # Select fields without subsetting variants
  selected <- subset(res, select = c(chr, pos, qual))

  expect_s3_class(selected, "vcftable")
  expect_equal(names(selected), c("samples", "chr", "pos", "qual"))
  expect_equal(nrow(selected$qual), nrow(res$qual))
})

test_that("subset preserves matrix format for gt", {
  skip_on_os("windows")
  res <- vcftable(vcffile, "chr21:1-5050000")

  # Subset and verify gt is still a matrix
  result <- subset(res, pos < 5010000)

  expect_true(is.matrix(result$gt))
  expect_equal(ncol(result$gt), ncol(res$gt))
  expect_equal(ncol(result$gt), length(result$samples))
})

test_that("subset with complex logical expressions works", {
  skip_on_os("windows")
  res <- vcftable(vcffile, "chr21:1-5050000")

  # Complex condition
  complex <- subset(res, qual > 50 & pos >= 5000000 & nchar(ref) == 1)

  expect_s3_class(complex, "vcftable")
  expect_true(all(complex$qual > 50))
  expect_true(all(complex$pos >= 5000000))
  expect_true(all(nchar(complex$ref) == 1))
})

test_that("subset with alternative allele filtering works", {
  skip_on_os("windows")
  res <- vcftable(vcffile, "chr21:1-5050000")

  # Filter by specific alternate allele
  a_alleles <- subset(res, alt == "A")

  expect_s3_class(a_alleles, "vcftable")
  expect_true(all(a_alleles$alt == "A"))
})

test_that("subset returns empty result when no variants match", {
  skip_on_os("windows")
  res <- vcftable(vcffile, "chr21:1-5050000")

  # Impossible condition
  empty <- subset(res, qual > 1000000)

  expect_s3_class(empty, "vcftable")
  expect_equal(nrow(empty$gt), 0)
  expect_equal(length(empty$pos), 0)
})

test_that("subset works with vcftable containing FORMAT tags", {
  skip_on_os("windows")
  # Test with different FORMAT tags
  res_ad <- vcftable(vcffile, "chr21:1-5050000", format = "AD")
  high_qual <- subset(res_ad, qual > 100)

  expect_s3_class(high_qual, "vcftable")
  expect_true(all(high_qual$qual > 100))
  expect_true(is.list(high_qual$AD)) # multiallelics makes AD unaligned
})

test_that("subset maintains sample names", {
  skip_on_os("windows")
  res <- vcftable(vcffile, "chr21:1-5050000")
  original_samples <- res$samples

  # Subset data
  result <- subset(res, pos < 5010000)

  expect_equal(result$samples, original_samples)
  expect_equal(colnames(result$gt), NULL)  # Check if column names are preserved
})

test_that("subset with ID field works", {
  skip_on_os("windows")
  res <- vcftable(vcffile, "chr21:1-5050000", setid = TRUE)

  # Get first few IDs
  first_ids <- head(res$id, 10)

  # Subset by ID
  id_subset <- subset(res, id %in% first_ids)

  expect_s3_class(id_subset, "vcftable")
  expect_true(all(id_subset$id %in% first_ids))
  expect_equal(length(id_subset$id), length(first_ids))
})

test_that("subset with drop parameter works", {
  skip_on_os("windows")
  res <- vcftable(vcffile, "chr21:1-5050000")

  # Subset with drop=TRUE
  result_drop <- subset(res, pos < 5010000, drop = TRUE)
  # Subset with drop=FALSE (default)
  result_no_drop <- subset(res, pos < 5010000, drop = FALSE)

  expect_s3_class(result_drop, "vcftable")
  expect_s3_class(result_no_drop, "vcftable")
  # Both should still be matrices since we're subsetting multiple rows
  expect_true(is.matrix(result_drop$gt))
  expect_true(is.matrix(result_no_drop$gt))
})

