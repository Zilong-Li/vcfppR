#' @title
#' Make sensible and beautiful plots based on various objects in vcfppR
#' 
#' @param obj object returned by vcfcomp 
#'
#' @param what what statisitcs to be plotted
#'
#' @param which.sample which sample among all to be plotted
#'
#' @param variant which types of variant are desired
#'
#' @param pop file contains population information
#'
#' @param ... parameters passed to graphics
#'
#' @export
vcfplot <- function(obj,
                    what = "r2",
                    which.sample = 1,
                    variant = c("SNP","INDEL"),
                    pop = NULL,
                    ...) {
  if(!is(obj, "vcfcomp") & !is(obj, "vcfsummary"))
    stop("the input object is not vcfcomp or vcfsummary class")
  colorpalette("colorblind")
  
  if(is(obj, "vcfcomp")){
    obj.names <- names(obj)
    if("r2" %in% obj.names & what == "r2")
      plot_r2(obj$r2, ...)
    if("pse" %in% obj.names & what == "pse")
      plot_pse(obj$pse[which.sample], ...)
    if("f1" %in% obj.names & what == "f1")
      plot_mat(obj$f1, which.sample, ...)
    if("nrc" %in% obj.names & what == "nrc")
      plot_mat(obj$nrc, which.sample, ...)
  }

  if(is(obj, "vcfsummary")){
    if(is.null(pop)) {
      svs <- obj$summary[-1]
      barplot(svs, ylim = c(0, 1.1*max(svs)),...)
    } else {
      # get labels and do plottiing
      ped <- read.table(pop, header=TRUE)
      i <- grep("Super|Population", colnames(ped))
      if(length(i)>1) 
        i <- grep("Super", colnames(ped))[1]
      upops <- unique(ped[,i])
      mat <- do.call(rbind, obj[-c(1,2)])
      out <- sapply(upops, function(pop) {
        id <- ped[ped[,i]==pop, ]$Sample
        ord <- match(id, obj$samples)
        mat[variant,ord]
      })
      boxplot(out, ...)
    }
  }  
}


plot_mat <- function(d, which.sample,
                     rm.na = TRUE, bin = NULL,
                     ylim = c(0,1), main = "",
                     xlab = "Minor allele frequency %",
                     ylab = expression(italic(r^2)),
                     ...) {
  d <- do.call(rbind, d[,"concordance"])
  if(rm.na) d <- d[complete.cases(d[,which.sample]),]
  bins <- get_bin(d)
  x <- log10(bins)
  labels <- 100 * bins
  if(is.null(bin)) bin <- seq_along(x)
  plot(0, 0,
       col = "transparent",
       axes = FALSE,
       xlim = range(x[bin]),
       ylim = ylim,
       xlab = xlab,
       ylab = ylab,
       main = main)
  for(i in seq_along(which.sample))
    points(x=x[bin], y=d[bin, which.sample[i]], col = i, ...)
  axis(side = 1, at = x[bin], labels = labels[bin], cex.axis=2) 
  axis(side = 2, at = seq(0, 1, 0.2), cex.axis=2)
}


plot_r2 <- function(d, rm.na = TRUE, bin = NULL,
                    ylim = c(0,1), main = "",
                    xlab = "Minor allele frequency %",
                    ylab = expression(italic(r^2)),
                    ...) {
  if(rm.na) d <- d[complete.cases(d[,"concordance"]),]
  bins <- get_bin(d)
  x <- log10(bins)
  labels <- 100 * bins
  if(is.null(bin)) bin <- seq_along(x)
  plot(0, 0,
       col = "transparent",
       axes = FALSE,
       xlim = range(x[bin]),
       ylim = ylim,
       xlab = xlab,
       ylab = ylab,
       main = main)
  points(x=x[bin], y=d[bin, "concordance"], ...)
  axis(side = 1, at = x[bin], labels = labels[bin], cex.axis=2) 
  axis(side = 2, at = seq(0, 1, 0.2), cex.axis=2)
}


plot_pse <- function(pse, extra = 10, col = 2, xaxt = "s",
                     xlab = "Genomic coordinates", ylab = "",
                     main = "", at = NULL,...) {
  COL <- palette()[col]
  pos <- lapply(pse, function(p) split_coordinates(p$pos))
  xmax <- max(unlist(pos)) + extra
  xmin <- max(c(0, min(unlist(pos))))
  plot(0, 0, col = "white", axes=FALSE,
       xlim = c(xmin, xmax), ylim = c(0, length(pse)),
       xlab = xlab, ylab = ylab, main = main)
  if(is.null(at)) at <- pretty(c(xmin, xmax))
  if(xaxt!="n")
    axis(1, at = at, ...)
  ## axis(1, at = axTicks(1))
  for(n in seq_along(pos)) {
    a <- c(xmin, pos[[n]], xmax) ## pad xmin and xmax
    for(l in 2:length(a)) {
      rect(a[l-1], n-1, a[l], n, col = ifelse(l %% 2 == 0, COL, add_alpha(COL, 0.5)))
    }
  }
}

