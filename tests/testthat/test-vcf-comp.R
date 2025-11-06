library(testthat)

rawvcf <- system.file("extdata", "raw.gt.vcf.gz", package="vcfppR")
imputedvcf <- system.file("extdata", "imputed.gt.vcf.gz", package="vcfppR")

svtruthvcf <- system.file("extdata", "platinum.sv.vcf.gz", package="vcfppR")
svuppvcf <- system.file("extdata", "svupp.call.vcf.gz", package="vcfppR")


test_that("can work with gtgq (GQ stratified concordance)", {
  skip_on_os(c("windows"), arch = NULL)
  truth <- vcftable(svtruthvcf)
  expect_equal(length(truth$id),1793L)
  truth$neighbors <-as.integer(sub(".*NumNeighbors=([^;]+).*", "\\1", truth$info))
  truth <- subset(truth, neighbors == 0) ## subset biallelic SVs
  expect_equal(length(truth$id),1048L)
  res <- vcfcomp(svuppvcf, truth, stats = "gtgq", region = "chr1")
  expect_equal(nrow(res$gtgq),7336L)
  expect_true(max(res$gtgq$meandisc) < 0.1)
})

test_that("can work for correlation r2 between DS and GT", {
  skip_on_os(c("windows"), arch = NULL)
  samples <- "HG00673,NA10840"
  res <- vcfcomp(imputedvcf, imputedvcf, bins = c(0,1),
                 by.sample = TRUE, samples = samples, stats = "all",
                 setid = TRUE)
  expect_equal(as.numeric(unlist(res$r2[1,])), c(15, 6, 1, 1), tolerance=1e-6)
  expect_identical(unlist(res$f1[[3]]), rep(1,2))
  expect_identical(unlist(res$nrc[[3]]), rep(1,2))
})

test_that("can work for r2 with af", {
  skip_on_os(c("windows"), arch = NULL)
  samples <- "HG00673,NA10840"
  d1 <- vcftable(imputedvcf,setid = T, info=F)
  af <- runif(15)
  names(af) <- d1$id
  affile <- tempfile()
  saveRDS(af, affile)
  res <- vcfcomp(imputedvcf, rawvcf, stats = "r2", bins = c(0,1),
                 af = affile, samples = samples, setid = TRUE)
  expect_identical(res[[2]][,3], 1)
})

test_that("can work for F1 score", {
  skip_on_os(c("windows"), arch = NULL)
  samples <- "HG00673,NA10840"
  res <- vcfcomp(rawvcf, imputedvcf, stats = "f1", samples = samples,
                 by.sample = TRUE, bins = c(0,1), setid = TRUE)
  expect_identical(paste0(res$samples, collapse = ","), samples)  
  expect_identical(unlist(res[[2]][[3]]), rep(1,2))
})

test_that("can work for NRC rate by sample", {
  skip_on_os(c("windows"), arch = NULL)
  samples <- "HG00673,NA10840"
  res <- vcfcomp(rawvcf, imputedvcf, stats = "nrc", samples = samples,
                 by.sample = TRUE, bins = c(0,1), setid = TRUE)
  expect_identical(paste0(res$samples, collapse = ","), samples)  
  expect_identical(unlist(res[[2]][[3]]), rep(1,2))
})

test_that("can work for PSE", {
  skip_on_os(c("windows"), arch = NULL)
  samples <- "HG00151,HG00380"
  (res <- vcfcomp(imputedvcf, imputedvcf, stats = "pse", samples = samples, setid = TRUE))
  expect_true(is.nan(res[[2]][[1]]$pse))
  expect_true((res[[2]][[2]]$pse==0))
  expect_true((res[[2]][[2]]$disc==0))
  expect_true(is.null(res[[2]][[2]]$pos))
  expect_identical(paste0(res$samples, collapse = ","), samples)
})



test_that("can work for vcftable objects", {
  skip_on_os(c("windows"), arch = NULL)
  samples <- "HG00673,NA10840"
  test <- vcftable(rawvcf, samples = samples, setid = T)
  truth <- vcftable(imputedvcf, samples = samples, setid = T)
  res <- vcfcomp(test, truth, stats = "nrc", samples = samples,
                 by.sample = TRUE, bins = c(0,1), setid = TRUE)
  expect_identical(paste0(res$samples, collapse = ","), samples)
  expect_identical(unlist(res[[2]][[3]]), rep(1,2))
})

test_that("can work with flip parameter", {
  skip_on_os(c("windows"), arch = NULL)
  samples <- "HG00673,NA10840"
  # Without flip
  res_noflip <- vcfcomp(imputedvcf, rawvcf, stats = "r2",
                        samples = samples, bins = c(0,1),
                        setid = TRUE, flip = FALSE)
  # With flip
  res_flip <- vcfcomp(imputedvcf, rawvcf, stats = "r2",
                      samples = samples, bins = c(0,1),
                      setid = TRUE, flip = TRUE)

  expect_s3_class(res_noflip, "vcfcomp")
  expect_s3_class(res_flip, "vcfcomp")
  # Results should be different when flip is used
  # (unless the data happens to be symmetric)
})

test_that("can work with by.variant parameter", {
  skip("to be fixed")
  ## skip_on_os(c("windows"), arch = NULL)
  samples <- "HG00673,NA10840"
  # by.variant = TRUE
  res <- vcfcomp(imputedvcf, rawvcf, stats = "r2",
                 samples = samples, bins = c(0,1),
                 by.variant = TRUE, setid = TRUE)

  expect_s3_class(res, "vcfcomp")
  expect_true("r2" %in% names(res))
  # Check structure is different from by.sample
  expect_true(is.matrix(res$r2) | is.data.frame(res$r2))
})

test_that("can work with AF file in text format", {
  skip_on_os(c("windows"), arch = NULL)
  samples <- "HG00673,NA10840"
  d1 <- vcftable(imputedvcf, setid = TRUE, info = FALSE)

  # Create text AF file
  af_data <- data.frame(
    chr = d1$chr,
    pos = d1$pos,
    ref = d1$ref,
    alt = d1$alt,
    af = runif(length(d1$chr))
  )
  affile <- tempfile()
  write.table(af_data, affile, row.names = FALSE, quote = FALSE)

  # Test with text AF file
  res <- vcfcomp(imputedvcf, rawvcf, stats = "r2", bins = c(0,1),
                 af = affile, samples = samples, setid = TRUE)

  expect_s3_class(res, "vcfcomp")
  expect_identical(res[[2]][,3], 1)

  unlink(affile)
})

test_that("can work with out parameter to save intermediate files", {
  skip_on_os(c("windows"), arch = NULL)
  samples <- "HG00673,NA10840"

  # Create temporary directory for output
  outdir <- tempdir()
  outprefix <- file.path(outdir, "test_vcfcomp")

  res <- vcfcomp(imputedvcf, rawvcf, stats = "r2",
                 samples = samples, bins = c(0,1),
                 setid = TRUE, out = outprefix)

  # Check that output files were created
  expect_true(file.exists(paste0(outprefix, ".af.rds")))
  expect_true(file.exists(paste0(outprefix, ".test.rds")))
  expect_true(file.exists(paste0(outprefix, ".truth.rds")))

  # Verify files can be read
  af_saved <- readRDS(paste0(outprefix, ".af.rds"))
  expect_true(is.numeric(af_saved) | is.null(af_saved))

  # Clean up
  unlink(paste0(outprefix, ".af.rds"))
  unlink(paste0(outprefix, ".test.rds"))
  unlink(paste0(outprefix, ".truth.rds"))
})

test_that("can work with return_pse_sites parameter", {
  skip_on_os(c("windows"), arch = NULL)
  samples <- "HG00151,HG00380"

  # With return_pse_sites = TRUE
  res_with_sites <- vcfcomp(imputedvcf, imputedvcf, stats = "pse",
                            samples = samples, setid = TRUE,
                            return_pse_sites = TRUE)

  # With return_pse_sites = FALSE (default)
  res_no_sites <- vcfcomp(imputedvcf, imputedvcf, stats = "pse",
                          samples = samples, setid = TRUE,
                          return_pse_sites = FALSE)

  expect_s3_class(res_with_sites, "vcfcomp")
  expect_s3_class(res_no_sites, "vcfcomp")

  # Check PSE components exist
  expect_true("pse" %in% names(res_with_sites))
  expect_true("pse" %in% names(res_no_sites))
})

test_that("can work with choose_random_start for PSE", {
  skip_on_os(c("windows"), arch = NULL)
  samples <- "HG00151,HG00380"

  # Test with choose_random_start = TRUE
  res_random <- vcfcomp(imputedvcf, imputedvcf, stats = "pse",
                        samples = samples, setid = TRUE,
                        choose_random_start = TRUE)

  # Test with choose_random_start = FALSE (default)
  res_fixed <- vcfcomp(imputedvcf, imputedvcf, stats = "pse",
                       samples = samples, setid = TRUE,
                       choose_random_start = FALSE)

  expect_s3_class(res_random, "vcfcomp")
  expect_s3_class(res_fixed, "vcfcomp")
  expect_true("pse" %in% names(res_random))
  expect_true("pse" %in% names(res_fixed))
})

test_that("can work with names parameter for sample renaming", {
  skip("to be fixed" )
  samples <- "HG00673,NA10840"
  new_names <- c("Sample1", "Sample2")

  res <- vcfcomp(imputedvcf, rawvcf, stats = "r2",
                 samples = samples, bins = c(0,1),
                 names = new_names, setid = TRUE)

  expect_s3_class(res, "vcfcomp")
  expect_identical(res$samples, new_names)
})

test_that("can work with stats='all'", {
  skip_on_os(c("windows"), arch = NULL)
  samples <- "HG00673,NA10840"

  res <- vcfcomp(imputedvcf, imputedvcf, bins = c(0,1),
                 by.sample = TRUE, samples = samples, stats = "all",
                 setid = TRUE)

  expect_s3_class(res, "vcfcomp")
  expect_true("r2" %in% names(res))
  expect_true("f1" %in% names(res))
  expect_true("nrc" %in% names(res))
  expect_identical(unlist(res$f1[[3]]), rep(1,2))
  expect_identical(unlist(res$nrc[[3]]), rep(1,2))
})

test_that("can work with both by.sample and by.variant FALSE", {
  skip_on_os(c("windows"), arch = NULL)
  samples <- "HG00673,NA10840"

  # Neither by.sample nor by.variant
  res <- vcfcomp(imputedvcf, rawvcf, stats = "r2",
                 samples = samples, bins = c(0,1),
                 by.sample = FALSE, by.variant = FALSE,
                 setid = TRUE)

  expect_s3_class(res, "vcfcomp")
  expect_true("r2" %in% names(res))
})

test_that("can handle different format combinations", {
  skip_on_os(c("windows"), arch = NULL)
  samples <- "HG00673,NA10840"

  # DS vs DS
  res_ds_ds <- vcfcomp(imputedvcf, imputedvcf,
                       formats = c("DS", "DS"),
                       stats = "r2", samples = samples,
                       bins = c(0,1), setid = TRUE)

  # GT vs GT
  res_gt_gt <- vcfcomp(rawvcf, rawvcf,
                       formats = c("GT", "GT"),
                       stats = "r2", samples = samples,
                       bins = c(0,1), setid = TRUE)

  expect_s3_class(res_ds_ds, "vcfcomp")
  expect_s3_class(res_gt_gt, "vcfcomp")
})
