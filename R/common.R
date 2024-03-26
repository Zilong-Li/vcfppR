
concordance_by_freq <- function(truthG, testDS, breaks, af, FUN,
                                which_snps = NULL,
                                flip = FALSE,
                                per_snp = FALSE,
                                per_ind = FALSE) {
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
  x <- cut(af, breaks = breaks, include.lowest = TRUE)
  if (per_ind) {
    cors_per_af <- tapply(1:length(x), x, function(w) {
      list(
        n = length(w),
        nA = sum(truthG[w, ], na.rm = TRUE),
        concordance = unlist(sapply(1:ncol(truthG), function(ind) {
          FUN(truthG[w, ind], testDS[w, ind])
        }))
      )
    })
  } else if (ncol(truthG) > 1 && per_snp) {
    # for multiple sample, calculate r2 per snp then average them
    cors_per_af <- tapply(1:length(x), x, function(w) {
      c(
        n = length(w),
        nA = sum(truthG[w, ], na.rm = TRUE),
        concordance = mean(sapply(w, function(ww) {
          FUN(truthG[ww, ], testDS[ww, ])
        }), na.rm = TRUE)
      )
    })
  } else {
    cors_per_af <- tapply(1:length(x), x, function(w) {
      c(
        n = length(w),
        nA = sum(truthG[w, ], na.rm = TRUE),
        concordance = FUN(truthG[w,], testDS[w,])
      )
    })
  }
  # fill with NA for AF bins without SNPs
  cors_per_af <- t(sapply(cors_per_af, function(a) {
    if (is.null(a[1])) {
      return(c(n = NA, nA = NA, concordance = NA))
    }
    a
  }))
  return(cors_per_af)
}

## sugar r2
R2 <- function(a, b) {
  cor(as.vector(a), as.vector(b), use = "pairwise.complete")**2
}

## follow hap.py
## truth\imputed
##        0      1      2
## 0      ignore FP     FP
## 1      FN     TP     FP/FN
## 2      FN     FP/FN     TP
## f1 = 2 * TP / (2 * TP + FP + FN)
F1 <- function(a, b) {
  o <- table(as.vector(a), as.vector(b), useNA = "always")
  ## make table square
  if(nrow(o)!=ncol(o)){
    if(nrow(o) == ncol(o)+1){
      o <- o[-nrow(o),]
    } else if (nrow(o)+1==ncol(o)){
      o <- o[,-ncol(o)]
    } else{
      warning("ONLY homozygous (0) found in either truth or test data")
      return(NA)
    }
  }
  if(all(dim(o)==c(4,4)))
    o <- o[1:3,1:3]
  if(all(dim(o)!=c(3,3))) {
    warning("F1 should be used only for a sample with genotypes of all types, hom ref(0), het(1) and hom alt(2)")
    return(NA)
  }
  TP <- o[2,2] + o[3,3]
  FP <- o[1,2] + o[1, 3] + o[2, 3] + o[3,2]
  FN <- o[2,1] +o[2,3] + o[3,1] + o[3,2]
  res <- 2 * TP / (2 * TP + FP + FN)
  return(res)
}

## follow GLIMPSE2_concordance
## None Reference Concordnace = 1 - (e0 + e1 + e2) / (e0 + e1 + e2 + m1 + m2)
## truth\imputed
##        0      1      2
## 0      ignore e0    e0
## 1      e1     m1    e1
## 2      e2     e2  m2
## a <- c(1, 2, 0, 1,1)
## b <- c(1, 1, 0, 0,1)
## NRC(a, b)
## 
NRC <- function(a, b) {
  o <- table(as.vector(a), as.vector(b), useNA = "always")
  ## make table square
  if(nrow(o)!=ncol(o)){
    if(nrow(o) == ncol(o)+1){
      o <- o[-nrow(o),]
    } else if (nrow(o)+1==ncol(o)){
      o <- o[,-ncol(o)]
    } else{
      warning("ONLY homozygous (0) found in either truth or test data")
      return(NA)
    }
  }
  if(all(dim(o)==c(4,4)))
    o <- o[1:3,1:3]
  if(all(dim(o)!=c(3,3))) {
    warning("NRC should be used only for a sample with genotypes of all types, hom ref(0), het(1) and hom alt(2)")
    return(NA)
  }
  mismatches <-  sum(c(o[1,2:3], o[2,1], o[2,3], o[3,1:2]))
  matches <- sum(c(o[2,2], o[3,3]))
  res <- mismatches / (mismatches+matches)
  return(1-res)
}

