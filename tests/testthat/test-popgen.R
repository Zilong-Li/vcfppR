vcffile <- system.file("extdata", "raw.gt.vcf.gz", package="vcfppR")

test_that("both vcftable and vcfreader work for calculating heterzygosity", {
  vcf <- vcftable(vcffile, vartype = "snps")
  res1 <- as.integer(colSums(vcf[["gt"]]==1, na.rm = TRUE))
  br <- vcfreader$new(vcffile)
  res2 <- rep(0L, br$nsamples())
  while(br$variant()) {
    if(br$isSNP()) {
      gt <- br$genotypes(TRUE) == 1
      gt[is.na(gt)] <- FALSE
      res2 <- res2 + gt
    }
  }
  expect_identical(as.integer(res1), res2)
})

test_that("heterzygosity only counts snps", {
  vcf <- vcftable(vcffile, vartype = "snps", pass = TRUE)
  res1 <- as.integer(colSums(vcf[["gt"]]==1, na.rm = TRUE))
  vcf <- vcfpopgen(vcffile, pass = TRUE, fun = "heterozygosity")
  res2 <- vcf$hets
  expect_identical(res1, res2)
})
