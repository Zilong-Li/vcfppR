#' @title
#' Compare two VCF/BCF files reporting various statistics
#'
#' @details
#' \code{vcfcomp} implements various statistics to compare two VCF/BCF files,
#' e.g. report genotype concordance, correlation stratified by allele frequency.
#'
#' @param test path to the comparison file (test), which can be a VCF/BCF file, vcftable object or saved RDS file.
#'
#' @param truth path to the baseline file (truth), which can be a VCF/BCF file, vcftable object or saved RDS file.
#'
#' @param formats character vector. the FORMAT tags to extract for the test and truth respectively.
#'               default c("DS", "GT") extracts 'DS' of the test and 'GT' of the truth.
#'
#' @param stats character. the statistics to be calculated. Supports the following options:
#'              \describe{
#'                \item{"r2"}{the Pearson correlation coefficient squared (default)}
#'                \item{"f1"}{the F1-score, good balance between sensitivity and precision}
#'                \item{"nrc"}{the Non-Reference Concordance rate}
#'                \item{"pse"}{the Phasing Switch Error rate}
#'                \item{"all"}{calculate r2, f1, and nrc together}
#'                \item{"gtgq"}{genotype quality-based concordance analysis}
#'                \item{"gtdp"}{depth-based concordance analysis}
#'              }
#' 
#' @param by.sample logical. calculate sample-wise concordance, which can be stratified by MAF bin.
#' 
#' @param by.variant logical. calculate variant-wise concordance, which can be stratified by MAF bin.
#'                  If both by.sample and by.variant are FALSE, then do calculations for all samples and variants together in a bin.
#' 
#' @param flip logical. flip the ref and alt variants
#' 
#' @param names character vector. reset samples' names in the test VCF.
#' 
#' @param bins numeric vector. break statistics into allele frequency bins.
#'              If NULL (default), bins are automatically generated with fine resolution for rare variants
#'              and coarser resolution for common variants (ranging from 0 to 0.5).
#'
#' @param af file path with allele frequency or a RDS file with a saved object for af.
#'           Format of the text file: a space-separated text file with five columns and a header named 'chr' 'pos' 'ref' 'alt' 'af'.
#'           If NULL, allele frequencies are calculated from the truth genotypes.
#'
#' @param out output prefix for saving objects into RDS file. If provided, creates three files:
#'            <out>.af.rds, <out>.test.rds, and <out>.truth.rds
#'
#' @param choose_random_start logical. choose random start for stats="pse". Defaults to FALSE.
#'
#' @param return_pse_sites logical. return phasing switch error sites when stats="pse". Defaults to FALSE.
#'
#' @param ... additional options passed to \code{vcftable}, such as 'samples', 'region', or 'pass'.
#'
#' @return a list object of class "vcfcomp" containing:
#'         \describe{
#'           \item{samples}{character vector of sample names}
#'           \item{[stats]}{the calculated statistics, named according to the 'stats' parameter.
#'                          For stats="all", returns r2, f1, and nrc components.}
#'         }
#' 
#' @author Zilong Li \email{zilong.dk@gmail.com}
#'
#' @examples
#' library('vcfppR')
#' # site-wise comparision stratified by allele frequency
#' test <- system.file("extdata", "imputed.gt.vcf.gz", package="vcfppR")
#' truth <- system.file("extdata", "raw.gt.vcf.gz", package="vcfppR")
#' samples <- "HG00673,NA10840"
#' res <- vcfcomp(test, truth, stats="r2", bins=c(0,1), samples=samples, setid=TRUE)
#' str(res)
#' 
#' # sample-wise comparision stratified by sample-level metrice e.g GQ
#' test <- system.file("extdata", "svupp.call.vcf.gz", package="vcfppR")
#' truth <- system.file("extdata", "platinum.sv.vcf.gz", package="vcfppR")
#' res <- vcfcomp(test, truth, stats = "gtgq", region = "chr1")
#' str(res)
#' 
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
  if( stats %in% c("r2", "nrc", "f1", "pse", "all")) {
    comp_stats_vs_maf(test, truth, formats, stats, by.sample, by.variant, flip, names, bins, af, out, choose_random_start, return_pse_sites, ...)
  } else {
    comp_stats_per_sample(test, truth, stats, ...)    
  }
}

## +--------------------+
## | internal functions |
## +--------------------+

## per site statisitcs
comp_stats_vs_maf <- function(test, truth,
                              formats,
                              stats,
                              by.sample,
                              by.variant,
                              flip,
                              names,
                              bins,
                              af,
                              out,
                              choose_random_start,
                              return_pse_sites,
                              ...) {
  
  if((stats=="f1" | stats == "nrc" | stats == "pse") & (formats[1] != "GT") & (stats != "all")) {
    message("stats F1 or NRC or PSE only uses GT format")
    formats[1] <- "GT"
  }

  if(is.null(bins)){
    bins <- sort(unique(c(
      c(0, 0.01 / 1000, 0.02 / 1000, 0.05 / 1000),
      c(0, 0.01 / 100, 0.02 / 100, 0.05 / 100),
      c(0, 0.01 / 10, 0.02 / 10, 0.05 / 10),
      c(0, 0.01 / 1, 0.02 / 1, 0.05 / 1),
      seq(0.1, 0.5, length.out = 5)
    )))
  }

  collapse <- ifelse(stats=="pse", FALSE, TRUE)

  if(is(truth, "vcftable")) {
    d1 <- test
  } else {
    d1 <- tryCatch( { suppressWarnings(readRDS(test)) }, error = function(e) {
      vcftable(test, format = formats[1], collapse = collapse, ...)
    } )
  }

  if(is(truth, "vcftable")) {
    d2 <- truth
  } else {
    d2 <- tryCatch( { suppressWarnings(readRDS(truth)) }, error = function(e) {
      vcftable(truth, format = formats[2], collapse = collapse, ...)
    } )
  }

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
    ## R2
    res.r2 <- concordance_by_freq(gt, ds, bins, af, R2, which_snps = sites,
                                  flip = flip, per_ind = by.sample, per_snp = by.variant)
    if(stats == "r2")
      return(list(samples = d1$samples, r2=res.r2))
    ## F1 and NRC
    d1 <- vcftable(test, format = "GT", ...)
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

## per sample statisitcs
comp_stats_per_sample <- function(test, truth, stats, ...) {
  formats <- switch(stats,
                    gtgq = c("GT", "GQ"),
                    gtdp = c("GT", "DP"),
                    stop("invalid option of stats"))

  if(is(test, "vcftable")) {
    d1 <- test
  } else {
    d1 <- tryCatch( { suppressWarnings(readRDS(test)) }, error = function(e) {
      d <- vcftable(test, format = formats[1], ...)
      d[[formats[2]]] <- vcftable(test, format = formats[2], ...)[[10]]
      d
    } )
  }

  if(is(truth, "vcftable")) {
    d2 <- truth
  } else {
    d2 <- tryCatch( { suppressWarnings(readRDS(truth)) }, error = function(e) {
      vcftable(truth, format = formats[1], ...)
    } )
  }

  w1 <- match(d2$id, d1$id)
  w2 <- match(d2$samples, d1$samples)
  if(sum(!d1$id[w1] == d2$id)>0) stop('Variant IDs do not match')
  if(sum(!d1$samples[w2] == d2$samples)>0) stop('Sample IDs do not match')
  gt <- d1[[10]][w1,w2]
  gq <- d1[[11]][w1,w2]
  call <- as.data.frame(cbind(
    'gq'=as.vector(gq),
    'disc'= as.vector(gt != d2[[10]])
  ))
  w <- !is.na(call$disc)
  call <- call[w,] ## remove NAs
  # then order by gq
  call <- call[order(call$gq, decreasing = T),]
  # return ranked mean discordance
  call$meandisc <- cumsum(call$disc)/1:sum(w)
  ret <- list(samples = d2$samples, stats = call)
  names(ret) <- c("samples", stats)
  class(ret) <- "vcfcomp"
  return(ret)
}
