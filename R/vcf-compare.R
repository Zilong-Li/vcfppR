#' @title
#' Compare two VCF/BCF files reporting various statistics
#'
#' @details
#' \code{vcfcomp} implements various statisitcs to compare two VCF/BCF files,
#' e.g. report genotype concocrdance, correlation stratified by allele frequency.
#' 
#' @param test path to the first VCF/BCF file referred as test.
#'
#' @param truth path to the second VCF/BCF file referred as truth.
#'
#' @param region region to subset in bcftools-like style: "chr1", "chr1:1-10000000"
#'
#' @param samples samples to subset in bcftools-like style.
#'                comma separated list of samples to include (or exclude with "^" prefix).
#'                e.g. "id01,id02", "^id01,id02".
#'
#' @param names rename samples in test VCF.
#' 
#' @param format character vector. the FORMAT tag to extract for comparison.
#'               default c("DS", "GT") is used to extract DS of target and GT of truth respectively.
#'
#' @param stats choose the statistics to be returned. eg. "r2", "f1"
#'
#' @param bins numeric vector. allele frequency bins to stratify with.
#'
#' @param af file path to allele frequency.
#'
#' @param out output prefix for saving objects into RDS file
#'
#' @param vartype restrict to specific type of variants. supports "snps","indels", "sv", "multisnps","multiallelics"
#' @param ids  character vector. restrict to sites with ID in the given vector.
#' @param qual logical. restrict to variants with QUAL > qual.
#'
#' @param pass logical. restrict to variants with FILTER = "PASS".
#'
#' @return return various statistics
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
vcfcomp <- function(test, truth, region = "", samples = "-", names = NULL,
                    format = c("DS", "GT"), stats = "r2", bins = NULL, af = NULL, out = NULL,
                    vartype = "snps", ids = NULL, qual = 0, pass = FALSE) {
  if(is.null(bins)){
    bins <- sort(unique(c(
      c(0, 0.01 / 1000, 0.02 / 1000, 0.05 / 1000),
      c(0, 0.01 / 100, 0.02 / 100, 0.05 / 100),
      c(0, 0.01 / 10, 0.02 / 10, 0.05 / 10),
      c(0, 0.01 / 1, 0.02 / 1, 0.05 / 1),
      seq(0.1, 0.5, length.out = 5)
    )))
  }
  if(stats=="f1" & format[1] != "GT") stop("F1 score is using GT format. please use format=c('GT','GT')")
  d1 <- vcftable(test, region = region, samples = samples, setid = TRUE, info = FALSE, format = format[1],
                 vartype = vartype, ids = NULL, qual = qual, pass = pass)
  if(!is.null(names) & is.vector(names)) d1$samples <- names
  samples <- paste0(d1$samples, collapse = ",") 
  d2 <- tryCatch( { suppressWarnings(readRDS(truth)) }, error = function(e) {
    vcftable(test, region = region, samples = samples, setid = TRUE, info = FALSE, format = format[2],
             vartype = vartype, ids = NULL, qual = qual, pass = pass)
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
  res <- NULL
  if(stats=="r2") {
    if(is.null(af)){
      af <- rowMeans(gt, na.rm = TRUE) / 2
    } else {
      af <- af[match(sites, af[,"id"]), "af"]
    }
    names(af) <- sites
    res <- r2_by_freq(bins, af, gt, ds, which_snps = sites, flip = FALSE)
    return(list(samples = d1$samples, r2=res))
  }
  if(stats=="f1"){
    res <- F1(gt, ds)
    return(list(samples = d1$samples, f1=res))
  }
}

