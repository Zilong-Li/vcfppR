test_that("clearInfo removes all INFO tags", {
  library(vcfppR)
  vcffile <- system.file("extdata", "raw.gt.vcf.gz", package="vcfppR")
  outvcf <- tempfile(fileext = ".vcf.gz")

  # Read VCF and write to output with INFO cleared
  br <- vcfreader$new(vcffile, "chr21:1-5050000")
  br$output(outvcf)

  # Process first variant
  expect_true(br$variant())

  # Check that INFO column has content before clearing
  info_before <- br$info()
  expect_true(nchar(info_before) > 0)

  # Clear INFO
  br$clearInfo()

  # Check that INFO column is now empty (should be ".")
  info_after <- br$info()
  expect_equal(info_after, "")

  # Write the modified variant
  br$write()

  # Process second variant without clearing INFO
  expect_true(br$variant())
  info_second <- br$info()
  expect_true(nchar(info_second) > 0)
  br$write()

  br$close()

  # Now verify the output file
  br2 <- vcfreader$new(outvcf)

  # First variant should have empty INFO
  expect_true(br2$variant())
  expect_equal(br2$info(), "")

  # Second variant should still have INFO
  expect_true(br2$variant())
  expect_true(nchar(br2$info()) > 0)

  # Clean up
  unlink(outvcf)
})

test_that("clearInfo works with multiple variants", {
  library(vcfppR)
  vcffile <- system.file("extdata", "raw.gt.vcf.gz", package="vcfppR")
  outvcf <- tempfile(fileext = ".vcf.gz")

  # Read VCF and clear INFO for all variants
  br <- vcfreader$new(vcffile, "chr21:1-5050000")
  br$output(outvcf)

  variant_count <- 0
  while(br$variant()) {
    br$clearInfo()
    br$write()
    variant_count <- variant_count + 1
  }
  br$close()

  expect_true(variant_count > 0)

  # Verify all variants have empty INFO
  br2 <- vcfreader$new(outvcf)
  while(br2$variant()) {
    expect_equal(br2$info(), "")
  }

  # Clean up
  unlink(outvcf)
})

test_that("clearInfo handles variants with no INFO", {
  library(vcfppR)
  # Create a VCF with minimal INFO
  outvcf <- tempfile(fileext = ".vcf.gz")
  bw <- vcfwriter$new(outvcf, "VCF4.1")
  bw$addContig("chr21")
  bw$addFORMAT("GT", "1", "String", "Genotype")
  bw$addSample("sample1")

  # Write variant with empty INFO (.)
  variant_line <- "chr21\t100\t.\tA\tT\t100\tPASS\t.\tGT\t0/1"
  bw$writeline(variant_line)
  bw$close()

  # Read and try to clear INFO (should not error)
  br <- vcfreader$new(outvcf)
  outvcf2 <- tempfile(fileext = ".vcf.gz")
  br$output(outvcf2)

  expect_true(br$variant())
  expect_equal(br$info(), "")

  # Should not error when clearing already empty INFO
  expect_silent(br$clearInfo())
  br$write()
  br$close()

  # Clean up
  unlink(outvcf)
  unlink(outvcf2)
})
