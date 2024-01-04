library(testthat)

vcffile <- system.file("extdata", "raw.gt.vcf.gz", package="vcfppR")
svfile <- system.file("extdata", "sv.vcf.gz", package="vcfppR")

test_that("vcfreader: constructor with vcfffile only", {
  br <- vcfreader$new(vcffile)
  expect_true(br$variant())
  expect_identical(br$chr(), "chr21")
  expect_identical(br$pos(), 5030082L)
  expect_identical(br$id(), "chr21:5030082:G:A")
  expect_identical(br$ref(), "G")
  expect_identical(br$alt(), "A")
})

test_that("vcfreader: constructor with vcfffile and region", {
  br <- vcfreader$new(vcffile, "chr21:5030082-")
  expect_true(br$variant())
  expect_identical(br$chr(), "chr21")
  expect_identical(br$pos(), 5030082L)
  expect_identical(br$id(), "chr21:5030082:G:A")
  expect_identical(br$ref(), "G")
  expect_identical(br$alt(), "A")
})

test_that("vcfreader: constructor with vcfffile, region and samples", {
  br <- vcfreader$new(vcffile, "chr21:5030082-", "HG00097,HG00099")
  expect_identical(br$samples(), c("HG00097", "HG00099"))
  expect_true(br$variant())
  expect_identical(br$genotypes(F), rep(0L, 4))
  expect_identical(br$chr(), "chr21")
  expect_identical(br$pos(), 5030082L)
  expect_identical(br$id(), "chr21:5030082:G:A")
  expect_identical(br$ref(), "G")
  expect_identical(br$alt(), "A")
})

test_that("vcfreader: get FORMAT item with Float type", {
  br <- vcfreader$new(vcffile, "chr21:5030347-", "HG00097,HG00099")
  expect_identical(br$samples(), c("HG00097", "HG00099"))
  expect_true(br$variant())
  expect_identical(br$pos(), 5030347L)
  expect_identical(all(is.na(br$formatFloat("AB"))), TRUE)
})

test_that("vcfreader: get FORMAT item with Integer type", {
  br <- vcfreader$new(vcffile, "chr21:5030347-", "HG00097,HG00099")
  expect_identical(br$samples(), c("HG00097", "HG00099"))
  expect_true(br$variant())
  expect_identical(br$pos(), 5030347L)
  expect_identical(br$formatInt("AD"), c(4L, 0L, 16L, 0L))
})

test_that("vcfreader: get FORMAT item with String type", {
  br <- vcfreader$new(svfile, "chr21:5114000-", "HG00096,HG00097")
  expect_identical(br$samples(), c("HG00096", "HG00097"))
  expect_true(br$variant())
  expect_identical(br$pos(), 5114000L)
  expect_identical(br$formatStr("EV"), c("RD","RD"))
})

test_that("vcfreader: get FORMAT item for the unexpected and error", {
  br <- vcfreader$new(vcffile, "chr21:5030347-", "HG00097,HG00099")
  br$variant()
  expect_identical(br$pos(), 5030347L)
  expect_error(br$formatInt("AA")) ## AA not exists
  expect_identical(br$formatInt("AD"), c(4L, 0L, 16L, 0L))
  ## use formatFloat will return unexpected values if AD in FORMAT if Integer
  expect_error(expect_identical(br$formatFloat("AD"), c(4L, 0L, 16L, 0L))) ## AA exists but of Integer type
})

test_that("vcfreader: get variants type", {
  br <- vcfreader$new(vcffile)
  i <- 0
  s <- 0
  m <- 0
  ms <- 0
  while(br$variant()) {
    if(br$isIndel()) i <- i + 1
    if(br$isSV()) s <- s + 1
    if(br$isMultiAllelics()) m <- m + 1
    if(br$isMultiAllelicSNP()) ms <- ms + 1
  }
  expect_identical(i, 6)
  expect_identical(s, 0)
  expect_identical(m, 4)
  expect_identical(ms, 3)
})

test_that("vcfreader: test variants type", {
  br <- vcfreader$new(vcffile)
  i1 <- 0
  i2 <- 0
  i3 <- 0
  i4 <- 0
  i5 <- 0
  i6 <- 0
  i7 <- 0
  i8 <- 0
  while(br$variant()) {
    if(br$hasSNP()) i1 <- i1 + 1
    if(br$hasINDEL()) i2 <- i2 + 1
    if(br$hasINS()) i3 <- i3 + 1
    if(br$hasDEL()) i4 <- i4 + 1
    if(br$hasMNP()) i5 <- i5 + 1
    if(br$hasBND()) i6 <- i6 + 1
    if(br$hasOTHER()) i7 <- i7 + 1
    if(br$hasOVERLAP()) i8 <- i8 + 1
  }
  expect_identical(i1, 66)
  expect_identical(i2, 6)
  expect_identical(i3, 3)
  expect_identical(i4, 4)
  expect_identical(i5, 0)
  expect_identical(i6, 0)
  expect_identical(i7, 0)
  expect_identical(i8, 1)
})

test_that("vcfreader: reading variant only", {
  br <- vcfreader$new(vcffile)
  br$variant()
  expect_identical(br$chr(), "chr21")
  expect_identical(br$pos(), 5030082L)
  expect_identical(br$id(), "chr21:5030082:G:A")
  expect_identical(br$ref(), "G")
  expect_identical(br$alt(), "A")
  expect_equal(br$qual(), 70.06, tolerance = 1e-6)
  expect_identical(br$filter(), "VQSRTrancheSNP99.80to100.00")
  expect_identical(br$info(), "AC=2;AF=0.000616523;AN=3244;DP=2498;FS=0;MLEAC=1;MLEAF=0.0003083;MQ=17.07;MQ0=0;QD=17.52;SOR=3.258;VQSLOD=-32.6;culprit=MQ;AN_EUR=658;AN_EAS=548;AN_AMR=520;AN_SAS=624;AN_AFR=894;AF_EUR=0;AF_EAS=0;AF_AMR=0.00384615;AF_SAS=0;AF_AFR=0;AC_EUR=0;AC_EAS=0;AC_AMR=2;AC_SAS=0;AC_AFR=0;AC_Het_EUR=0;AC_Het_EAS=0;AC_Het_AMR=0;AC_Het_SAS=0;AC_Het_AFR=0;AC_Het=0;AC_Hom_EUR=0;AC_Hom_EAS=0;AC_Hom_AMR=2;AC_Hom_SAS=0;AC_Hom_AFR=0;AC_Hom=2;HWE_EUR=1;ExcHet_EUR=1;HWE_EAS=1;ExcHet_EAS=1;HWE_AMR=0.00192678;ExcHet_AMR=1;HWE_SAS=1;ExcHet_SAS=1;HWE_AFR=1;ExcHet_AFR=1;HWE=0.000308356;ExcHet=1;ME=0;AN_EUR_unrel=508;AN_EAS_unrel=448;AN_AMR_unrel=330;AN_SAS_unrel=484;AN_AFR_unrel=666;AF_EUR_unrel=0;AF_EAS_unrel=0;AF_AMR_unrel=0;AF_SAS_unrel=0;AF_AFR_unrel=0;AC_EUR_unrel=0;AC_EAS_unrel=0;AC_AMR_unrel=0;AC_SAS_unrel=0;AC_AFR_unrel=0;AC_Het_EUR_unrel=0;AC_Het_EAS_unrel=0;AC_Het_AMR_unrel=0;AC_Het_SAS_unrel=0;AC_Het_AFR_unrel=0;AC_Hom_EUR_unrel=0;AC_Hom_EAS_unrel=0;AC_Hom_AMR_unrel=0;AC_Hom_SAS_unrel=0;AC_Hom_AFR_unrel=0")
  br$rmInfoTag("AC")
  print(br$info())
  print(br$infoInt("AC")) ## TODO: AC value still exists
  af <- br$infoFloat("AF")
  ## expect_equal(0.00061652297, 0.00061652300, tolerance = 1e-6)
  expect_equal(af, 0.000616523, tolerance = 1e-6)
  dp <- br$infoInt("DP")
  expect_identical(dp, 2498L)
  gt <- br$genotypes(TRUE)
  expect_identical(length(gt), 3202L)
  gt <- br$genotypes(FALSE)
  expect_identical(length(gt), 6404L)
  pl <- br$formatInt("PL")
  expect_identical(length(pl), 9606L)
})

test_that("vcfreader: read and output variant", {
  br <- vcfreader$new(vcffile)
  br$variant()
  br$setCHR("chr22")
  expect_identical(br$chr(), "chr22")
  br$setPOS(9999L)
  expect_identical(br$pos(), 9999L)
  br$setID("rs122334")
  expect_identical(br$id(), "rs122334")
  expect_identical(br$ref(), "G")
  expect_identical(br$alt(), "A")
  br$setRefAlt("A,GA")
  expect_identical(br$ref(), "A")
  expect_identical(br$alt(), "GA")
  expect_identical(br$infoInt("AC"), 2L)
  br$setInfoInt("AC", 3L)
  expect_identical(br$infoInt("AC"), 3L)
  br$setInfoFloat("AF", 0.3)
  expect_equal(br$infoFloat("AF"), 0.3, tolerance  = 1e-6)
  br$setInfoStr("VariantType", "indel")
  expect_identical(br$infoStr("VariantType" ), "indel")
  ## output current variant to another vcf
  s1 <- br$string()
  outvcf <- paste0(tempfile(), ".vcf.gz")
  br$output(outvcf)
  br$write()
  br$close()
  br <- vcfreader$new(outvcf)
  br$variant()
  s2 <- br$string()
  expect_identical(s1, s2)
})


