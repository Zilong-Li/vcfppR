
r2_by_freq <- function(breaks, af, truthG, testDS, which_snps = NULL, flip = FALSE, per_snp = FALSE, per_ind = FALSE) {
  if (!is.null(which_snps)) {
    af <- af[which_snps]
    truthG <- truthG[which_snps, ]
    testDS <- testDS[which_snps, ]
  }
  truthG <- as.matrix(truthG)
  testDS <- as.matrix(testDS)
  af <- as.numeric(af)
  if (flip) {
    w <- af > 0.5
    af[w] <- 1 - af[w]
    truthG[w, ] <- 2 - truthG[w, ]
    testDS[w, ] <- 2 - testDS[w, ]
  }
  x <- cut(af, breaks = breaks)
  if (per_ind) {
    cors_per_af <- tapply(1:length(x), x, function(w) {
      c(
        n = length(w),
        nA = sum(truthG[w, ], na.rm = TRUE),
        simple = mean(sapply(1:ncol(truthG), function(ind) {
          cor(truthG[w, ind], testDS[w, ind], use = "pairwise.complete")**2
        }), na.rm = TRUE),
        norm = NA
      )
    })
  } else if (ncol(truthG) > 1 && per_snp) {
    # for multiple sample, calculate r2 per snp then average them
    cors_per_af <- tapply(1:length(x), x, function(w) {
      c(
        n = length(w),
        nA = sum(truthG[w, ], na.rm = TRUE),
        simple = mean(sapply(w, function(ww) {
          cor(truthG[ww, ], testDS[ww, ], use = "pairwise.complete")**2
        }), na.rm = TRUE),
        norm = mean(sapply(w, function(ww) {
          cor(truthG[ww, ] - 2 * af[ww], testDS[ww, ] - 2 * af[ww], use = "pairwise.complete")**2
        }), na.rm = TRUE)
      )
    })
  } else {
    cors_per_af <- tapply(1:length(x), x, function(w) {
      c(
        n = length(w),
        nA = sum(truthG[w, ], na.rm = TRUE),
        simple = cor(as.vector(truthG[w, ]), as.vector(testDS[w, ]), use = "pairwise.complete")**2,
        norm = cor(as.vector(truthG[w, ] - 2 * af[w]), as.vector(testDS[w, ] - 2 * af[w]), use = "pairwise.complete")**2
      )
    })
  }
  # fill with NA for AF bins without SNPs
  cors_per_af <- t(sapply(cors_per_af, function(a) {
    if (is.null(a[1])) {
      return(c(n = NA, nA = NA, simple = NA, norm = NA))
    }
    a
  }))
  return(cors_per_af)
}


F1 <- function(a, b) {
  stopifnot(dim(a)==dim(b))
  sapply(seq_len(ncol(a)), function(i) {
    o <- table(a[,i], b[,i])
    if(sum(colnames(o) == c("0", "1", "2")) != 3) {
      warning("F1 should be used only for a sample with genotypes of all types, hom ref(0), het(1) and hom alt(2)")
      return(NA)
    }
    TP <- o[2,2] + o[3,3]
    FP <- o[1,2] + o[1, 3] + o[2, 3] + o[3,2]
    FN <- o[2,1] +o[2,3] + o[3,1] + o[3,2]
    return(2 * TP / (2 * TP + FP + FN))
  })
}

