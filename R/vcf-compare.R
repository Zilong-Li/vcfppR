#' @title
#' Compare two VCF/BCF files reporting various statistics
#'
#' @details
#' \code{vcfcomp} implements various statisitcs to compare two VCF/BCF files,
#' e.g. report genotype concocrdance, correlation stratified by allele frequency.
#' 
#' @param test path to the first VCF/BCF file referred as test.
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
#' 
#' @param bysample logical. calculate concordance for each samples, then average by bins.
#' 
#' @param byvariant logical. calculate concordance for each variant, then average by bins.
#'                  if both bysample and by variant are TRUE, then do average on all samples first.
#'                  if both bysample and by variant are FALSE, then do average on all samples and variants.
#' 
#' @param names character vector. reset samples' names in test VCF.
#' 
#' @param bins numeric vector. break statistics into allele frequency bins.
#'
#' @param af file path to allele frequency text file or saved RDS file.
#'
#' @param out output prefix for saving objects into RDS file
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
#' res <- vcfcomp(test, truth, stats="f1", format=c('GT','GT'), samples=samples)
#' str(res)
#' @export
vcfcomp <- function(test, truth,
                    formats = c("DS", "GT"),
                    stats = "r2",
                    bysample = FALSE,
                    byvariant = FALSE,
                    names = NULL,
                    bins = NULL,
                    af = NULL,
                    out = NULL,
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
  if((stats=="f1" | stats == "nrc") & formats[1] != "GT") {
    message("F1 score or NRC rate use GT format")
    formats[1] <- "GT"
  }
  d1 <- vcftable(test, format = formats[1], setid = TRUE, ...)
  if(!is.null(names) & is.vector(names)) d1$samples <- names
  samples <- paste0(d1$samples, collapse = ",") 
  d2 <- tryCatch( { suppressWarnings(readRDS(truth)) }, error = function(e) {
    vcftable(test, format = formats[2], setid = TRUE, ...)
  } )
  sites <- intersect(d1$id,  d2$id)
  ## chr pos ref alt af
  if(!is.null(af)){
    af <- tryCatch( { suppressWarnings(readRDS(af)) }, error = function(e) {
      af <- read.table(af, header = TRUE)
      af$id <- paste0(af[,"chr"], "_", af[,"pos"], "_", af[,"ref"], "_", af[,"alt"])
      subset(af, select = c(id, af))
    } )
    sites <-  intersect(af[,"id"], sites) ## use intersect sites only
  }
  ## save some useful objects
  if(!is.null(out)){
    saveRDS(af, file.path(paste0(out, ".af.rds")))
    saveRDS(truth, file.path(paste0(out, ".truth.rds")))
  }
  ord <- match(d1$samples, d2$samples)
  ds <- d1[[10]]
  ds <- ds[match(sites, d1$id), ]
  gt <- d2[[10]]
  gt <- gt[match(sites, d2$id), ord]
  rownames(gt) <- sites
  rownames(ds) <- sites
  if(is.null(af)){
    af <- rowMeans(gt, na.rm = TRUE) / 2
  } else {
    af <- af[match(sites, af[,"id"]), "af"]
  }
  names(af) <- sites
  res <- NULL
  if(stats=="r2") {
    res <- concordance_by_freq(gt, ds, bins, af, R2, which_snps = sites,
                               flip = FALSE, per_ind = bysample, per_snp = byvariant)
    return(list(samples = d1$samples, r2=res))
  }
  if(stats=="f1"){
    res <- concordance_by_freq(gt, ds, bins, af, F1, which_snps = sites,
                               flip = FALSE, per_ind = bysample, per_snp = byvariant)
    return(list(samples = d1$samples, f1=res))
  }
  if(stats=="nrc"){
    res <- concordance_by_freq(gt, ds, bins, af, NRC, which_snps = sites,
                               flip = FALSE, per_ind = bysample, per_snp = byvariant)
    return(list(samples = d1$samples, nrc=res))
  }
}

