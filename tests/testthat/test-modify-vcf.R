library(testthat)

test_that("modify the genotypes", {
  outvcf <- paste0(tempfile(), ".vcf.gz")
  bw <- vcfwriter$new(outvcf, "VCF4.3")
  bw$addContig("chr20")
  bw$addINFO("AF", "A", "Float", "Estimated allele frequency in the range (0,1)");
  bw$addFORMAT("GT", "1", "String", "Genotype");
  bw$addSample("NA12878")
  bw$addSample("NA12879")
  s1 <- "chr20\t2006060\trs146931526\tG\tC\t100\tPASS\tAF=0.000998403\tGT\t1|0\t1/1"
  bw$writeline(s1)
  bw$close()
  ## tests
  br <- vcfreader$new(outvcf)
  br$variant() ## get a variant record
  (g0 <- br$genotypes(F))
  ## if you wanna change the phasing of genotypes,
  ## call setPhasing before setGenotypes
  br$setPhasing(c(1L, 1L)) 
  g1 <- c(1L, 1L, NA, 0L)
  br$setGenotypes(g1)
  (g2 <- br$genotypes(F))
  expect_false(isTRUE(all.equal(g0, g2)))
  expect_identical(g1, g2)
})
