#' @title
#' Make sensible and beautiful plots based on various objects in vcfppR
#' 
#' @param obj object returned by vcftable, vcfcomp, vcfsummary
#'
#' @param which.sample which sample to be plotted. NULL will aggregate all samples.
#'
#' @param which.format which FORMAT field to be plotted. Defaults will use the 10-th names.
#'
#' @param variant which types of variant are desired
#'
#' @param pop file contains population information
#'
#' @param ... parameters passed to graphics
#'
#' @export
vcfplot <- function(obj,
                    which.sample = NULL,
                    which.format = 10,
                    variant = c("SNP","INDEL"),
                    pop = NULL,
                    ...) {
  if(!(is(obj, "vcfcomp") | is(obj, "vcfsummary") | is(obj, "vcftable")))
    stop("the input object is not vcftable, vcfcomp or vcfsummary class")
  colorpalette("colorblind")
  
  if(is(obj, "vcfcomp")){
    obj.names <- names(obj)
    if("pse" %in% obj.names)
      plot_pse(obj$pse[which.sample], ...)
    else if("r2" %in% obj.names)
      plot_comp(obj, "r2", which.sample, ...)
    else if("f1" %in% obj.names)
      plot_comp(obj, "f1", which.sample, ...)
    else if("nrc" %in% obj.names)
      plot_comp(obj, "nrc", which.sample, ...)
    else
      message("no available plots for unnamed vcfcomp objects")
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

  if(is(obj, "vcftable")){
    plot_scatter(obj, which.sample, which.format, ...)
  }
}


plot_comp <- function(obj, stats, which.sample,
                      rm.na = TRUE, bin = NULL,
                      ylim = c(0,1),
                      ylab = NULL,
                      main = NULL,
                      xlab = "Minor allele frequency %",
                      ...) {
  bysample <- TRUE
  d <- tryCatch({ do.call(rbind, obj[[stats]][,"accuracy"])},
                error = function(e){
                  FALSE
                })
  if(!is.matrix(d)){
    bysample <- FALSE
    d <- obj[[stats]]
  }
  if(bysample & is.null(which.sample)){
    message("which.sample is NULL. will set which.sample as 1")
    which.sample <- 1
  }
  if(is.null(which.sample)){
    if(rm.na)
      d <- d[complete.cases(d[,"accuracy"]),]
  } else {
    if(rm.na)
      d <- d[complete.cases(d[,which.sample]),]
  }
  bins <- get_bin(d)
  x <- log10(bins)
  labels <- 100 * bins
  if(is.null(bin)) bin <- seq_along(x)
  if(is.null(ylab) && (stats!="r2")) ylab <- toupper(stats)
  if(is.null(ylab) && (stats=="r2")) ylab <- expression(italic(r^2))
  plot(0, 0,
       col = "transparent",
       axes = FALSE,
       xlim = range(x[bin]),
       ylim = ylim,
       xlab = xlab,
       ylab = ylab,
       main = main)
  if(is.null(which.sample)){
    points(x=x[bin], y=d[bin, "accuracy"], ...)
  } else {
    for(i in seq_along(which.sample))
      points(x=x[bin], y=d[bin, which.sample[i]], col = i, ...)
  }
  axis(side = 1, at = x[bin], labels = labels[bin], cex.axis=2) 
  axis(side = 2, at = seq(0, 1, 0.2), cex.axis=2)
}


plot_pse <- function(pse, extra = 10, col = 2, xaxt = "s",
                     xlab = "Genomic coordinates", ylab = NULL,
                     main = NULL, at = NULL,...) {
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


plot_scatter <- function(obj,
                         which.sample,
                         which.format = 10,
                         xlab = "Genomic coordinates",
                         ylab = NULL,
                         main = NULL,
                         ...){
  stopifnot(is(obj, "vcftable"))
  if(!is.matrix(obj[[which.format]]))
    warning("the format is not well aligned across all samples")
  if(is.null(which.sample)){
    message("which.sample is NULL. will set which.sample as 1")
    which.sample <- 1
  }
  if(length(which.sample)>1){
    message("which.sample has length > 1. select the first one")
    which.sample <- which.sample[1]
  }
  x <- obj$pos
  y <- obj[[which.format]]
  objattr <- names(obj)
  if(is.null(ylab))
    ylab <- ""
  if(is.null(main))
    main <- objattr[which.format]
  plot(x, y[,which.sample],
       xlab = xlab,
       ylab = ylab,
       main = main,
       ...)
}
