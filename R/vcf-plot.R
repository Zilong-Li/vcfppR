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
    else if("gtgq" %in% obj.names)
      plot_call(obj, "gtgq", ...)
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

#' @title Plot variants on haplotypes across multiple samples
#'
#' @description
#' Visualizes variant positions and alleles on both haplotypes for multiple samples.
#' Each sample is represented by two horizontal tracks (one per haplotype), with variants
#' colored according to their type (SNP, insertion, deletion) and allele (reference or alternate).
#' Large gaps between variants can be automatically compressed for better visualization.
#'
#' @param vcffiles Character vector of VCF/BCF file paths or URLs. Each file represents one sample.
#' @param region Character string specifying the genomic region to visualize (e.g., "chr1:1000-5000").
#' @param types Character vector of variant types to include in the plot.
#'        Valid options are "SNP" (single nucleotide polymorphisms),
#'        "DEL" (deletions), and "INS" (insertions). Default: c("SNP", "DEL", "INS").
#' @param shrink_threshold Numeric value specifying the minimum gap size (in base pairs)
#'        between variants that will trigger compression. Gaps larger than this threshold
#'        are shrunk to improve visualization density. Default: 1000.
#' @param xlab Character string for the x-axis label. Default: "Genomic position".
#' @param ylab Character string for the y-axis label. Default: "Haplotypes of each sample".
#' @param main Character string for the plot title. Default: NULL (no title).
#' @param ... Additional graphical parameters passed to the base plot function.
#'
#' @details
#' The function reads variant data from multiple VCF files using \code{\link{vcftable}} with
#' \code{collapse=FALSE} to preserve haplotype phasing information. Each sample is displayed
#' as two horizontal tracks representing the two haplotypes (h1 and h2).
#'
#' Variant types are distinguished by color:
#' \itemize{
#'   \item SNPs: Green (reference allele) or Yellow (alternate allele)
#'   \item Deletions: Light blue (reference allele only)
#'   \item Insertions: Dark blue (reference allele) or Dark orange (alternate allele)
#' }
#'
#' When large gaps exist between variants (exceeding \code{shrink_threshold}), the function
#' compresses these regions and marks them with red dashed lines and "..." text to indicate
#' the compression. This feature helps visualize sparse variant distributions more effectively.
#'
#' @return Invisibly returns NULL. The function is called for its side effect of creating a plot.
#'
#' @seealso \code{\link{vcftable}}, \code{\link{vcfplot}}
#'
#' @examples
#' \dontrun{
#' # Plot variants from three samples in a specific region
#' vcf_files <- c("sample1.vcf.gz", "sample2.vcf.gz", "sample3.vcf.gz")
#' plot_variants_per_haplotype(vcf_files, region = "chr20:1000000-1100000")
#'
#' # Plot only SNPs and insertions with custom threshold
#' plot_variants_per_haplotype(vcf_files,
#'                            region = "chr20:1000000-1100000",
#'                            types = c("SNP", "INS"),
#'                            shrink_threshold = 5000)
#'
#' # Customize plot appearance
#' plot_variants_per_haplotype(vcf_files,
#'                            region = "chr20:1000000-1100000",
#'                            main = "Variant Distribution",
#'                            xlab = "Position (bp)",
#'                            cex.axis = 0.8)
#' }
#'
#' @export
plot_variants_per_haplotype <- function(vcffiles, region,
                                        types = c('SNP', 'DEL', 'INS'),
                                        shrink_threshold = 1000, 
                                        xlab = "Genomic position",
                                        ylab = "Haplotypes of each sample",
                                        main = NULL,
                                        ...) {
  dl <- lapply(vcffiles, function(f) {
    dat <- vcftable(f, region = region, collapse = F, info = F)
    d <- data.frame(pos = integer(), ref = character(), alt = character(), h1 = integer(), h2 = integer())
    # wif no vars
    if(length(dat$gt)>0) {
      gt <- dat$gt
      d <- data.frame(pos = dat$pos, ref = dat$ref, alt = dat$alt, h1 = gt[,1], h2 = gt[,2])
    }
    return(d)
  })

  # Get all variant positions across samples
  all_pos <- sort(unique(unlist(lapply(dl, function(d) d$pos))))
  
  if(length(all_pos) == 0) {
    # No variants to plot
    xlim <- c(0, 1000)
    plot(1, col = 'transparent', xlim = xlim, ylim = c(0, length(dl)+2), 
         axes = T, xlab = xlab, ylab = ylab, main = main, ...)
    return()
  }
  
  # Find large gaps between variants
  gaps <- diff(all_pos)
  large_gaps <- which(gaps > shrink_threshold)
  
  # Create mapping between original and compressed positions
  compressed_pos <- all_pos
  shrink_lines <- c()  # positions for vertical lines
  
  if(length(large_gaps) > 0) {
    cumulative_shrink <- 0
    
    for(i in large_gaps) {
      gap_size <- gaps[i]
      shrink_amount <- gap_size - shrink_threshold/2
      
      # Mark positions for vertical lines (before shrinking)
      line_pos1 <- all_pos[i] + shrink_threshold/4 - cumulative_shrink
      line_pos2 <- all_pos[i+1] - shrink_threshold/4 - cumulative_shrink - shrink_amount
      shrink_lines <- c(shrink_lines, line_pos1, line_pos2)
      
      # Adjust positions after this gap
      compressed_pos[(i+1):length(compressed_pos)] <- compressed_pos[(i+1):length(compressed_pos)] - shrink_amount
      cumulative_shrink <- cumulative_shrink + shrink_amount
    }
  }
  
  # Create position mapping function
  map_position <- function(pos) {
    if(length(pos) == 0) return(numeric(0))
    
    mapped <- pos
    if(length(large_gaps) > 0) {
      cumulative_shrink <- 0
      
      for(i in large_gaps) {
        gap_size <- gaps[i]
        shrink_amount <- gap_size - shrink_threshold/2
        gap_start <- all_pos[i]
        gap_end <- all_pos[i+1]
        
        # Positions after the gap start get shifted
        after_gap <- pos > gap_start
        mapped[after_gap] <- mapped[after_gap] - shrink_amount
        
        cumulative_shrink <- cumulative_shrink + shrink_amount
      }
    }
    return(mapped)
  }
  
  # Update variant positions in data
  dl_mapped <- lapply(dl, function(d) {
    if(nrow(d) > 0) {
      d$pos <- map_position(d$pos)
    }
    return(d)
  })
  
  xlim <- range(c(compressed_pos, shrink_lines)) + c(-10, 10)
  ylim <- c(0, length(dl)+1)

  colorpalette('colorblind')
  plot(1, col = 'transparent', xlim = xlim, ylim = ylim, axes = T, xlab = xlab, ylab = ylab, main = main, ...)
  
  # Add vertical lines to indicate shrinked regions (capped to data region)
  if(length(shrink_lines) > 0) {
    for(i in seq(1, length(shrink_lines), 2)) {
      # Draw capped vertical lines instead of full ablines
      segments(x0 = shrink_lines[i], y0 = 0, x1 = shrink_lines[i],
               y1 = length(dl), col = "red", lty = 2, lwd = 1)
      segments(x0 = shrink_lines[i+1], y0 = 0, x1 = shrink_lines[i+1],
               y1 = length(dl), col = "red", lty = 2, lwd = 1)

      # Add "..." text between lines
      mid_pos <- mean(c(shrink_lines[i], shrink_lines[i+1]))
      text(mid_pos, length(dl)/2, "...", cex = 1.2, col = "red", font = 2)
    }
  }
  
  # Build legend with proper color palette
  legend_colors <- c()
  legend_labels <- c()
  pal <- palette()  # Get current palette

  if('SNP' %in% types) {
    legend_colors <- c(legend_colors, pal[4], pal[5])
    legend_labels <- c(legend_labels, "SNP Ref", "SNP Alt")
  }
  if('DEL' %in% types) {
    legend_colors <- c(legend_colors, pal[3])
    legend_labels <- c(legend_labels, "DEL Ref")
  }
  if('INS' %in% types) {
    legend_colors <- c(legend_colors, pal[6], pal[7])
    legend_labels <- c(legend_labels, "INS Ref", "INS Alt")
  }

  if(length(shrink_lines) > 0) {
    legend_colors <- c(legend_colors, "red")
    legend_labels <- c(legend_labels, "Compressed")
  }

  # Improved legend with better layout
  ## n_items <- length(legend_labels)
  ## legend_ncol <- min(n_items, 4)  # Maximum 4 columns for readability

  legend("top", legend = legend_labels, fill = legend_colors, cex = 1.1,
         ncol = length(legend_labels), bty = "n", bg = "white", xpd = FALSE)

  for(iid in seq_along(dl_mapped)) {
    add_haplotypes_per_iid(dl_mapped[[iid]], iid-1, xlim[1], xlim[2], types)
  }
}


## +--------------------+
## | internal functions |
## +--------------------+

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

plot_call <- function(obj, stats,
                      xlab = NULL,
                      ylab = NULL,
                      main = NULL,
                      ...) {
  if(is.null(xlab)) xlab <- "# called genotypes"
  if(is.null(ylab)) ylab <- "Mean genotype discordance"
  if(is.null(main)) main <- paste0("stats='", stats,"'")
  ## plot(x, d[x,"meandisc"], bty = 'l', type = 'l')
  d <- obj[[stats]]

  if(stats == "gtgq") {
    changed <- which(diff(d[,1]) != 0)+1
    x <- changed[seq(1,length(changed), 4)] ## select every 4 
    plot(changed, d[changed,"meandisc"],
         xlab = xlab,
         ylab = ylab,
         main = main,
         ...)
  } else {
    stop("no plots for this. please submit a PR!")
  }
}



add_haplotypes_per_iid <- function(dd, iid, xleft, xright, types = c('SNP', 'DEL', 'INS')) {

  colorpalette('colorblind')
  
  refdel <- 3 # lightblue
  refins <- 6 # darkblue
  altsnp <- 5 # yellow
  refsnp <- 4 # green
  altins <- 7 # darkorange

  rect(xleft, iid, xright, iid+0.2, lty = 1, border = NA, col = 1)
  rect(xleft, iid+0.5, xright, iid+0.7, lty = 1, border = NA, col = 1)

  w <- nchar(dd$alt) < nchar(dd$ref) # DEL
  if(('DEL' %in% types) & (sum(w)>0)) {
    a <- dd[w,]
    ## s <- ifelse(a$h1 == 1, a$alt, a$ref)
    s <- a$ref
    bc <- rep('white', length(a$pos))
    bc[which(a$h1 == 0)] <- refdel
    rect(a$pos, iid, a$pos+nchar(s), iid+0.2, border = bc, col = bc)

    bc <- rep('white', length(a$pos))
    bc[which(a$h2 == 0)] <- refdel
    rect(a$pos, iid+0.5, a$pos+nchar(s), iid+0.7, border = bc, col = bc)
  }

  w <- nchar(dd$alt) > nchar(dd$ref) # INS
  if(('INS' %in% types) & (sum(w)>0)) {
    a <- dd[w,]

    ## s <- ifelse(a$h1 == 1, a$alt, a$ref)
    s <- a$alt
    bc <- rep(altins, length(a$pos))
    bc[which(a$h1 == 0)] <- refins
    rect(a$pos, iid, a$pos+1, iid+0.2, border = bc, col = bc)

    bc <- rep(altins, length(a$pos))
    bc[which(a$h1 == 0)] <- 'transparent'  # we don't plot ref
    rect(a$pos, iid+0.2, a$pos+nchar(s), iid+0.3, border = bc, col = bc)

    bc <- rep(altins, length(a$pos))
    bc[which(a$h2 == 0)] <- refins
    rect(a$pos, iid+0.5, a$pos+1, iid+0.7, border = bc, col = bc)

    bc <- rep(altins, length(a$pos))
    bc[which(a$h2 == 0)] <- 'transparent'  # we don't plot ref
    rect(a$pos, iid+0.7, a$pos+nchar(s), iid+0.8, border = bc, col = bc)
  }

  w <- nchar(dd$alt) == nchar(dd$ref) # SNP
  if(('SNP' %in% types) & (sum(w)>0)) {
    a <- dd[w,]
    bc <- rep(altsnp, length(a$pos))
    bc[which(a$h1 == 0)] <- refsnp
    rect(a$pos, iid, a$pos+1, iid+0.2, border = bc, col = bc)

    bc <- rep(altsnp, length(a$pos))
    bc[which(a$h2 == 0)] <- refsnp
    rect(a$pos, iid+0.5, a$pos+1, iid+0.7, border = bc, col = bc)
  }
}


