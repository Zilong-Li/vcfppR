#' @title
#' Compare two VCF/BCF files reporting various statistics
#'
#' @details
#' \code{vcfcomp} implements various statisitcs to compare two VCF/BCF files,
#' e.g. report genotype concocrdance, correlation stratified by allele frequency.
#' 
#' @param test path to the first VCF/BCF file referred as test, or saved RDS file.
#'
#' @param truth path to the second VCF/BCF file referred as truth, or saved RDS file.
#'
#' @param formats character vector. the FORMAT tags to extract for the test and truth respectively.
#'               default c("DS", "GT") extracts 'DS' of the target and 'GT' of the truth.
#'
#' @param stats the statistics to be calculated. supports the following.
#'              "r2": pearson correlation coefficient ** 2.
#'              "f1": F1-score, good balance between sensitivity and precision.
#'              "nrc": Non-Reference Concordance rate
#'              "pse": Phasing Switch Error rate
#' 
#' @param by.sample logical. calculate concordance for each samples, then average by bins.
#' 
#' @param by.variant logical. calculate concordance for each variant, then average by bins.
#'                  if both bysample and by variant are TRUE, then do average on all samples first.
#'                  if both bysample and by variant are FALSE, then do average on all samples and variants.
#' 
#' @param flip logical. flip the ref and alt variants
#' 
#' @param names character vector. reset samples' names in the test VCF.
#' 
#' @param bins numeric vector. break statistics into allele frequency bins.
#'
#' @param af file path to allele frequency text file or saved RDS file.
#' 
#' @param out output prefix for saving objects into RDS file
#'
#' @param choose_random_start choose random start for stats="pse"
#'
#' @param return_pse_sites boolean. return phasing switch error sites
#'
#' @param ... options passed to \code{vcftable}
#'
#' @return a list of various statistics
#' @author Zilong Li \email{zilong.dk@gmail.com}
#'
#' @examples
#' library('vcfppR')
#' test <- system.file("extdata", "imputed.gt.vcf.gz", package="vcfppR")
#' truth <- system.file("extdata", "imputed.gt.vcf.gz", package="vcfppR")
#' samples <- "HG00133,HG00143,HG00262"
#' res <- vcfcomp(test, truth, stats="f1", format=c('GT','GT'), samples=samples, setid=TRUE)
#' str(res)
#' @export
vcfcomp <- function(test, truth,
                    formats = c("DS", "GT"),
                    stats = "r2",
                    by.sample = FALSE,
                    by.variant = FALSE,
                    flip = FALSE,
                    names = NULL,
                    bins = NULL,
                    af = NULL,
                    out = NULL,
                    choose_random_start = FALSE,
                    return_pse_sites = FALSE,
                    ...) {
  if(is.null(bins)){
    bins <- sort(unique(c(
      c(0, 0.01 / 1000, 0.02 / 1000, 0.05 / 1000),
      c(0, 0.01 / 100, 0.02 / 100, 0.05 / 100),
      c(0, 0.01 / 10, 0.02 / 10, 0.05 / 10),
      c(0, 0.01 / 1, 0.02 / 1, 0.05 / 1),
      seq(0.1, 0.5, length.out = 5)
    )))
  }
  if((stats=="f1" | stats == "nrc" | stats == "pse") & (formats[1] != "GT") & (stats != "all")) {
    message("stats F1 or NRC or PSE only uses GT format")
    formats[1] <- "GT"
  }
  collapse <- ifelse(stats=="pse", FALSE, TRUE)
  d1 <- tryCatch( { suppressWarnings(readRDS(test)) }, error = function(e) {
    vcftable(test, format = formats[1], collapse = collapse, ...)
  } )
  d2 <- tryCatch( { suppressWarnings(readRDS(truth)) }, error = function(e) {
    vcftable(truth, format = formats[2], collapse = collapse, ...)
  } )
  if(!is.null(names) & is.vector(names)) d1$samples <- names
  ## chr pos ref alt af
  sites <- intersect(d1$id,  d2$id)
  if(!is.null(af)){
    af <- tryCatch( { suppressWarnings(readRDS(af)) }, error = function(e) {
      af <- read.table(af, header = TRUE)
      af$id <- paste0(af[,"chr"], "_", af[,"pos"], "_", af[,"ref"], "_", af[,"alt"])
      aaf <- af[,"af"]
      names(aaf) <- af[,"id"]
      aaf
    } )
    sites <-  intersect(names(af), sites) ## use intersect sites only
  }
  ## save some useful objects
  if(!is.null(out)){
    saveRDS(af, file.path(paste0(out, ".af.rds")))
    saveRDS(test, file.path(paste0(out, ".test.rds")))
    saveRDS(truth, file.path(paste0(out, ".truth.rds")))
  }
  ord <- match(d1$samples, d2$samples)
  if(is.na(sum(ord)))
    stop("the samples name in two VCF files is inconsistent. please set `names`")
  if(!collapse) ord <- c(sapply(ord, function(i) c(2*i-1, 2*i)))
  ds <- as.matrix(d1[[10]])
  ds <- as.matrix(ds[match(sites, d1$id), ])
  gt <- as.matrix(d2[[10]])
  gt <- as.matrix(gt[match(sites, d2$id), ord])
  rownames(gt) <- sites
  rownames(ds) <- sites
  if(is.null(af)){
    af <- rowMeans(gt, na.rm = TRUE) / 2
  } else {
    af <- af[match(sites, names(af))]
  }
  names(af) <- sites
  if(stats == "all") {
    ## F2
    res.r2 <- concordance_by_freq(gt, ds, bins, af, R2, which_snps = sites,
                                  flip = flip, per_ind = FALSE, per_snp = by.variant)
    if(stats == "r2")
      return(list(samples = d1$samples, r2=res.r2))
    ## F1 and NRC
    d1 <- vcftable(test, format = "GT", setid = TRUE, ...)
    ds <- d1[[10]]
    ds <- ds[match(sites, d1$id), ]
    rownames(ds) <- sites
    res.f1 <- concordance_by_freq(gt, ds, bins, af, F1, which_snps = sites,
                                  flip = flip, per_ind = by.sample, per_snp = by.variant)
    res.nrc <- concordance_by_freq(gt, ds, bins, af, NRC, which_snps = sites,
                                   flip = flip, per_ind = by.sample, per_snp = by.variant)
    ret <- list(samples = d1$samples, r2=res.r2, f1=res.f1, nrc=res.nrc)
  } else {
    res <- switch(stats,
                  pse = PSE(ds, gt, sites, choose_random_start, return_pse_sites),
                  r2 = concordance_by_freq(gt, ds, bins, af, R2, which_snps = sites,
                                           flip = flip, per_ind = by.sample, per_snp = by.variant),
                  f1 = concordance_by_freq(gt, ds, bins, af, F1, which_snps = sites,
                                           flip = flip, per_ind = by.sample, per_snp = by.variant),
                  nrc = concordance_by_freq(gt, ds, bins, af, NRC, which_snps = sites,
                                            flip = flip, per_ind = by.sample, per_snp = by.variant))
    ret <- list(d1$samples, res)
    names(ret) <- c("samples", stats)
  }
  class(ret) <- "vcfcomp"
  return(ret)
}

