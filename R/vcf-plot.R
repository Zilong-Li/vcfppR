#' @title
#' Do plotting based on various objects in vcfppR
#' 
#' @param obj object returned by vcfcomp 
#'
#' @export
vcfplot <- function(obj, which.sample = 1, ...) {
  stopifnot(is(r, "vcfcomp"))
  colorpalette("colorblind")
  obj.names <- names(obj)
  if("r2" %in% obj.names) plot_r2(obj$r2, ...)
  if("pse" %in% obj.names) plot_pse(obj$pse[which.sample])
  if("f1" %in% obj.names) plot_mat(obj$f1, which.sample, ...)
  if("nrc" %in% obj.names) plot_mat(obj$nrc, which.sample, ...)
}


plot_mat <- function(d, which.sample,rm.na = TRUE, w = NULL,
                     ylim = c(0,1), main = "",
                     xlab = "Minor allele frequency %",
                     ylab = expression(italic(r^2)),
                     ...) {
  d <- do.call(rbind, d[,"concordance"])
  if(rm.na) d <- d[complete.cases(d[,which.sample]),]
  bins <- get_bin(d)
  x <- log10(bins)
  labels <- 100 * bins
  if(is.null(w)) w <- seq_along(x)
  plot(0, 0,
       col = "transparent",
       axes = FALSE,
       xlim = range(x[w]),
       ylim = ylim,
       xlab = xlab,
       ylab = ylab,
       main = main)
  for(i in seq_along(which.sample))
    points(x=x[w], y=d[w, which.sample[i]], col = i, ...)
  axis(side = 1, at = x[w], labels = labels[w], cex.axis=2) 
  axis(side = 2, at = seq(0, 1, 0.2), cex.axis=2)
}


plot_r2 <- function(d, rm.na = TRUE, w = NULL, ylim = c(0,1),
                    main = "",
                    xlab = "Minor allele frequency %",
                    ylab = expression(italic(r^2)),
                    ...) {
  if(rm.na) d <- d[!sapply(d[, "concordance"], is.na), ]
  bins <- get_bin(d)
  x <- log10(bins)
  labels <- 100 * bins
  if(is.null(w)) w <- seq_along(x)
  plot(0, 0,
       col = "transparent",
       axes = FALSE,
       xlim = range(x[w]),
       ylim = ylim,
       xlab = xlab,
       ylab = ylab,
       main = main)
  points(x=x[w], y=d[w, "concordance"], ...)
  axis(side = 1, at = x[w], labels = labels[w], cex.axis=2) 
  axis(side = 2, at = seq(0, 1, 0.2), cex.axis=2)
}

