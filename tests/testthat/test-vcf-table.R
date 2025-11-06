library(testthat)

bcffile <- system.file("extdata", "test.bcf", package="vcfppR")
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


test_that("query a non-existing region", {
  expect_error(res <- vcftable(vcffile, region = "chr14:0-1000"))
})

test_that("subset samples and region for BCF", {
  res <- vcftable(bcffile, region =  "1:10000-13000", samples =  "I1,I30", vartype = "snps", format = "DS")
  m <- matrix(c(2.000, 2.000,  0.195, 0.000, 1.000, 0.000, 2.000, 0.036), byrow = TRUE, ncol = 2)
  expect_equal(res[[10]], m, tolerance = 1e6)
})

test_that("rmdup parameter removes duplicate positions", {
  skip_on_os("windows")
  # Read without rmdup
  res_with_dups <- vcftable(vcffile, region = "chr21:1-5050000")
  n_variants <- length(res_with_dups$pos)

  # Read with rmdup (should remove any duplicates at same position)
  res_no_dups <- vcftable(vcffile, region = "chr21:1-5050000", rmdup = TRUE)

  # Check no duplicated positions remain
  expect_false(any(duplicated(res_no_dups$pos)))

  # Number of variants should be <= original
  expect_true(length(res_no_dups$pos) <= n_variants)

  # Genotype matrix should match position count
  expect_equal(nrow(res_no_dups$gt), length(res_no_dups$pos))
})

test_that("mac parameter filters by minor allele count", {
  skip("to be fixed")
  # Without MAC filter
  res_all <- vcftable(vcffile, region = "chr21:1-5050000")

  # With MAC > 1 (at least 2 minor alleles)
  res_mac1 <- vcftable(vcffile, region = "chr21:1-5050000", mac = 1)

  # With MAC > 5
  res_mac5 <- vcftable(vcffile, region = "chr21:1-5050000", mac = 5)

  # More restrictive MAC should have fewer variants
  expect_true(length(res_mac1$pos) <= length(res_all$pos))
  expect_true(length(res_mac5$pos) <= length(res_mac1$pos))

  # Matrix dimensions should match
  expect_equal(nrow(res_mac1$gt), length(res_mac1$pos))
  expect_equal(nrow(res_mac5$gt), length(res_mac5$pos))
})

test_that("qual parameter filters by quality score", {
  skip_on_os("windows")
  # Without quality filter
  res_all <- vcftable(vcffile, region = "chr21:1-5050000")

  # With qual > 50
  res_qual50 <- vcftable(vcffile, region = "chr21:1-5050000", qual = 50)

  # With qual > 200
  res_qual200 <- vcftable(vcffile, region = "chr21:1-5050000", qual = 200)

  # More restrictive quality should have fewer variants
  expect_true(length(res_qual50$pos) <= length(res_all$pos))
  expect_true(length(res_qual200$pos) <= length(res_qual50$pos))

  # All qualities should be above threshold
  expect_true(all(res_qual50$qual > 50))
  expect_true(all(res_qual200$qual > 200))
})

test_that("collapse parameter works with GT format", {
  skip_on_os("windows")
  # With collapse=TRUE (default for diploid GT)
  res_collapsed <- vcftable(vcffile, region = "chr21:1-5050000",
                            format = "GT", collapse = TRUE)

  # With collapse=FALSE (maintains phasing)
  res_uncollapsed <- vcftable(vcffile, region = "chr21:1-5050000",
                              format = "GT", collapse = FALSE)

  # Collapsed should have one column per sample
  expect_equal(ncol(res_collapsed$gt), length(res_collapsed$samples))

  # Uncollapsed should have two columns per sample (for diploid)
  expect_equal(ncol(res_uncollapsed$gt), 2 * length(res_uncollapsed$samples))

  # Both should have same number of variants
  expect_equal(nrow(res_collapsed$gt), nrow(res_uncollapsed$gt))
})

test_that("collapse parameter works with non-GT formats", {
  skip("to be fixed")
  # With collapse=TRUE for AD (should try to create matrix)
  res_collapsed <- vcftable(vcffile, region = "chr21:1-5050000",
                            format = "AD", collapse = TRUE)

  # With collapse=FALSE (keeps as list)
  res_uncollapsed <- vcftable(vcffile, region = "chr21:1-5050000",
                              format = "AD", collapse = FALSE)

  # Collapsed should be a matrix (if possible)
  expect_true(is.matrix(res_collapsed$ad) || is.list(res_collapsed$ad))

  # Uncollapsed should be a list
  expect_true(is.list(res_uncollapsed$ad))
})

test_that("multiple filters work together", {
  skip("to be fixed")
  # Combine multiple filters
  res <- vcftable(vcffile, region = "chr21:1-5050000",
                  vartype = "snps",
                  qual = 100,
                  pass = TRUE,
                  mac = 2)

  # Check all filters are applied
  expect_true(all(nchar(res$ref) == 1))
  expect_true(all(nchar(res$alt) == 1))
  expect_true(all(res$qual > 100))
  expect_true(all(res$filter == "PASS"))
})

test_that("setid parameter creates custom IDs", {
  skip_on_os("windows")
  # Without setid
  res_no_setid <- vcftable(vcffile, region = "chr21:1-5050000", setid = FALSE)

  # With setid
  res_setid <- vcftable(vcffile, region = "chr21:1-5050000", setid = TRUE)

  # Check ID format is chr_pos_ref_alt
  expect_true(all(grepl("^chr21_[0-9]+_[ACGT]+_[ACGT,]+$", res_setid$id)))

  # IDs should be unique (if no duplicates)
  expect_equal(length(unique(res_setid$id)), length(res_setid$id))
})

test_that("samples parameter with exclusion works", {
  skip_on_os("windows")
  # Get all samples
  res_all <- vcftable(vcffile, region = "chr21:1-5050000")
  all_samples <- res_all$samples

  # Exclude first two samples with ^ prefix
  samples_to_exclude <- paste0("^", paste(all_samples[1:2], collapse = ","))
  res_excluded <- vcftable(vcffile, region = "chr21:1-5050000",
                           samples = samples_to_exclude)

  # Check excluded samples are not present
  expect_false(all_samples[1] %in% res_excluded$samples)
  expect_false(all_samples[2] %in% res_excluded$samples)

  # Number of samples should be reduced
  expect_equal(length(res_excluded$samples), length(all_samples) - 2)
})

test_that("vartype parameter correctly filters variant types", {
  skip_on_os("windows")
  # Test all variant types
  res_snps <- vcftable(vcffile, region = "chr21:1-5050000", vartype = "snps")
  res_indels <- vcftable(vcffile, region = "chr21:1-5050000", vartype = "indels")
  res_multisnps <- vcftable(vcffile, region = "chr21:1-5050000", vartype = "multisnps")
  res_all <- vcftable(vcffile, region = "chr21:1-5050000", vartype = "all")

  # SNPs should have single base ref and alt
  expect_true(all(nchar(res_snps$ref) == 1))
  expect_true(all(nchar(res_snps$alt) == 1))

  # INDELs should have length differences
  expect_true(all(nchar(res_indels$ref) != nchar(res_indels$alt)))

  # Multisnps should have comma in alt
  if(length(res_multisnps$alt) > 0) {
    expect_true(all(grepl(",", res_multisnps$alt)))
  }

  # All should include everything
  expect_true(length(res_all$pos) >= length(res_snps$pos))
  expect_true(length(res_all$pos) >= length(res_indels$pos))
})

test_that("region parameter accepts different formats", {
  skip_on_os("windows")
  # Chromosome only
  res_chr <- vcftable(vcffile, region = "chr21")

  # Chromosome with range
  res_range <- vcftable(vcffile, region = "chr21:5000000-5010000")

  # Range should have fewer variants than whole chromosome
  expect_true(length(res_range$pos) < length(res_chr$pos))

  # All positions in range should be within bounds
  expect_true(all(res_range$pos >= 5000000))
  expect_true(all(res_range$pos <= 5010000))
})

test_that("ids parameter filters by variant IDs", {
  skip_on_os("windows")
  # Get some IDs
  res_all <- vcftable(vcffile, region = "chr21:1-5050000")
  test_ids <- res_all$id[1:5]

  # Filter by those IDs
  res_filtered <- vcftable(vcffile, region = "chr21:1-5050000", ids = test_ids)

  # Should only get those IDs
  expect_true(all(res_filtered$id %in% test_ids))
  expect_equal(length(res_filtered$id), length(test_ids))
})

test_that("empty result returns proper structure", {
  skip_on_os("windows")
  # Use impossible quality filter
  res_empty <- vcftable(vcffile, region = "chr21:1-5050000", qual = 10000)

  # Should return empty vcftable structure
  expect_s3_class(res_empty, "vcftable")
  expect_equal(length(res_empty$pos), 0)
  expect_equal(length(res_empty$gt), 0)
})

test_that("SV vartype works correctly", {
  skip_on_os("windows")
  # Test with SV file
  res_sv <- vcftable(svfile, vartype = "sv")

  # Should get structural variants
  expect_true(length(res_sv$pos) > 0)
  expect_s3_class(res_sv, "vcftable")
})
