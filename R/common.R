
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

calculate_pse_2ways <- function(test,
                                truth,
                                which_snps,
                                i_option = 0,
                                seed = NULL) {
  ## for testing
  if (!is.null(seed)) set.seed(seed)
  truth <- truth[which_snps, ]
  test <- test[which_snps, ]
  which_sites <-
    rowSums(truth == 0 | truth == 1, na.rm = TRUE) == 2 &
    rowSums(truth, na.rm = TRUE) == 1 &
    rowSums(is.na(truth)) == 0
  truth <- matrix(truth[which_sites, ], ncol=2,byrow=T)
  test <- matrix(test[which_sites, ], ncol=2, byrow=T)
  snps <- which_snps[which_sites]
  if (nrow(test) == 0) {
    warning("No heterozygous sites!")
    return(NA)
  }
  ## as these sites, discrepency as well
  disc <- sum(rowSums(test) != 1)
  ## round test data for now. choose hets at random
  ## specifically remove from consideration double phase switch errors
  w <- rowSums(test) == 1
  w2 <- which(diff(abs(test[w, 1] - truth[w, 1])) != 0)
  to_remove <- NULL
  if (length(w2) > 0) {
    w3 <- which(diff(w2) == 1)
    if (length(w3) > 0) {
      for (a in w3) {
        c <- w2[c(a, a + 1)]
        to_remove <- c(to_remove, which(w)[c])
      }
    }
  }
  ## double pse are two consecutive
  if (length(to_remove) > 0) {
    test <- test[-to_remove, ]
    truth <- truth[-to_remove, ]
    snps <- snps[-to_remove]
  }
  ##
  if (i_option == 0) {
    ## only consider non-discrepent sites
    ## chose best start
    if (test[1, 1] != truth[1, 1]) {
      test <- test[, c(2, 1)]
    }
    ## calculate number of differences
    w <- rowSums(test) == 1
    if (sum(w) == 0) {
      warning("Test has no hets! possibly an error or homo over region")
      return(NA)
      ## switches1 <- cbind(i1 = NA, i2 = NA, l1 = NA, l2 = NA)
      phase_errors <- 0
      phase_sites <- 0
    } else {
      y <- diff(abs(test[w, 1] - truth[w, 1])) != 0
      snps <- snps[y]
      phase_errors <- sum(y)
      phase_sites <- sum(w) - 1
      ## we need rownames(test) <- 1:nrow(test) beforehand
      ## s <- as.integer(rownames(test[w, , drop = FALSE][c(as.logical(y), FALSE), , drop = FALSE]))
      ## e <- as.integer(rownames(test[w, , drop = FALSE][c(FALSE, as.logical(y)), , drop = FALSE]))
      ## switches1 <- cbind(i1 = s, i2 = e)
    }
  }
  if (i_option == 1) {
    choose_at_random <- which(rowSums(test) != 1)
    if (length(choose_at_random) > 0) {
      test[choose_at_random, ] <- 0
      r <- sample(
        c(1, 2),
        length(choose_at_random),
        replace = TRUE
      )
      test[cbind(choose_at_random, r)] <- 1
    }
    ## chose best start
    if (test[1, 1] != truth[1, 1]) {
      test <- test[, c(2, 1)]
    }
    ## calculate number of differences
    y <- diff(abs(test[, 1] - truth[, 1])) != 0
    snps <- snps[y]
    phase_errors <- sum(y)
    phase_sites <- nrow(test) - 1
  }
  ##
  return(
    list(
      values = c(
        phase_errors = phase_errors,
        phase_sites = phase_sites,
        disc_errors = disc,
        dist_n = nrow(test)
      ),
      sites = snps
    )
  )
}

PSE <- function(test, truth,
                sites = NULL,
                choose_random_start = FALSE,
                return_pse_sites = TRUE) {
  stopifnot(is.matrix(test) & is.matrix(truth))
  stopifnot(ncol(test) %% 2 == 0)
  N <- ncol(test) / 2
  if(is.null(sites)) sites <- 1:nrow(test)
  pos <- NULL
  o <- lapply(1:N, function(i) {
    itruth <-truth[,1:2+(i-1)*2] # truth
    itest <- test[,1:2+(i-1)*2] # test
    res <- calculate_pse_2ways(itest, itruth, which_snps = sites, i_option = choose_random_start)
    if(is.list(res)){
      values <- res$values
      pse <- values["phase_errors"] / values["phase_sites"]
      disc <- values["disc_errors"] / values["dist_n"]
      if(return_pse_sites) pos <- res$sites
      list(pse=pse, disc=disc, pos=pos)
    } else { 
      NA 
    }
  })
  return(o)
}
