library(testthat)

test_that("modify the genotypes", {
  ## skip_on_os(c("windows"), arch = NULL)
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
  br$close()
  ##  the original vcf can not be modified
  br <- vcfreader$new(outvcf)
  br$variant() ## get a variant record
  (g3 <- br$genotypes(F))
  expect_identical(g0, g3)
})

test_that("modify item in FORMAT", {
  ## skip_on_os(c("windows"), arch = NULL)
  ## creat a VCF with GP in FORMAT
  outvcf <- paste0(tempfile(), ".vcf.gz")
  bw <- vcfwriter$new(outvcf, "VCF4.3")
  bw$addContig("chr20")
  bw$addINFO("AF", "A", "Float", "Estimated allele frequency in the range (0,1)");
  bw$addFORMAT("GP", "3", "Float", "Posterior genotype probability of 0/0, 0/1, and 1/1");
  bw$addSample("NA12878")
  bw$addSample("NA12879")
  s1 <- "chr20\t2006060\trs146931526\tG\tC\t100\tPASS\tAF=0.000998403\tGP\t0.966,0.034,0\t0.003,0.872,0.125"
  bw$writeline(s1)
  bw$close()
  ## tests
  br <- vcfreader$new(outvcf)
  expect_true(br$variant()) ## get a variant record
  br$string()
  gp <- br$formatFloat("GP")
  gp <- array(gp, c(3, br$nsamples()))
  ds <- gp[2,] + gp[3,] * 2
  ## will prompt uses a message if `output` is not called 
  br$addFORMAT("DS", "1", "Float", "Diploid dosage")
  br$addINFO("INFO", "1", "Float", "INFO score of imputation")
  ## now open another file for output 
  newvcf <- paste0(tempfile(), ".vcf.gz")
  br$output(newvcf)
  ## add INFO, DS in header first
  br$addINFO("INFO", "1", "Float", "INFO score of imputation")
  br$addFORMAT("DS", "1", "Float", "Diploid dosage")
  br$addFORMAT("AC", "1", "Integer", "Allele counts")
  br$addFORMAT("STR", "1", "String", "Test String type")
  print(br$header())
  ## set DS in FORMAT now
  br$setFormatFloat("DS", ds)
  ## test if DS presents 
  expect_identical(br$formatFloat("DS"), ds)
  ## more tests
  br$setFormatInt("AC", c(3L, 4L))
  expect_false(br$setFormatStr("STR","HHH,JJJ")) ## length(s) %% nsamples != 0
  expect_true(br$setFormatStr("STR","HHHJJJ")) ## length(s) %% nsamples == 0
  print(br$string())
})
