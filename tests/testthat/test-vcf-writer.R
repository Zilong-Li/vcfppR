
test_that("vcfwriter: writing variant works", {
  outvcf <- paste0(tempfile(), ".vcf.gz")
  bw <- vcfwriter$new(outvcf, "VCF4.3")
  bw$addContig("chr20")
  bw$addINFO("AF", "A", "Float", "Estimated allele frequency in the range (0,1)");
  bw$addFORMAT("GT", "1", "String", "Genotype");
  bw$addSample("NA12878")
  s1 <- "chr20\t2006060\trs146931526\tG\tC\t100\tPASS\tAF=0.000998403\tGT\t1|0"
  bw$writeline(s1)
  bw$close()
  ## tests
  br <- vcfreader$new(outvcf)
  br$variant()
  s2 <- gsub("\n", "", br$string())
  expect_identical(br$chr(), "chr20")
  print(s1)
  print(s2)
  expect_identical(s1, s2)
})



