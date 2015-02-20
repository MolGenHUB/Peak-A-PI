# This should detect and install missing packages before loading them.
pkg <- c("shiny", "Rsamtools", "GenomicRanges", "GenomicAlignments")
npkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
if (length(npkg)) install.packages(npkg)

lapply(pkg, function(x) {
  library(x, character.only = T)
})

getSeqnamesBam <- function(file) {
  bf <- BamFile(file)
  seqnames(seqinfo(bf))
}

# getSeqlengthsBam <- function(file) { 
#   bf <- BamFile(file)
#   seqlengths(seqinfo(bf)) 
# }

getCoverage <- function(bam, chr) {
  if (!file.exists(paste(bam, "bai", sep = "."))) {
    stop("Unable to find index")
  }
  gal <- readGAlignments(bam, format = "BAM")
  gr <- sort(granges(gal[seqnames(gal) %in% chr]))
  rm(gal)
  
  len <- seqlengths(gr)[[chr]]
  ## seqinfo <- Seqinfo(chr,len, NA, NA)
  
  cvg.pg <- coverage(gr[strand(gr) == "+"])
  cvg.ng <- coverage(gr[strand(gr) == "-"])
  
  cvg.p5 <- coverage(resize(gr[strand(gr) == "+"], 1, "start"))
  cvg.p3 <- coverage(resize(gr[strand(gr) == "+"], 1, "end"))
  
  cvg.n5 <- coverage(resize(gr[strand(gr) == "-"], 1, "start"))
  cvg.n3 <- coverage(resize(gr[strand(gr) == "-"], 1, "end"))
  
  cvg.pg <- as.numeric(cvg.pg[[chr]])
  cvg.ng <- as.numeric(cvg.ng[[chr]])
  
  cvg.p5 <- as.numeric(cvg.p5[[chr]])
  cvg.p3 <- as.numeric(cvg.p3[[chr]])
  
  cvg.n5 <- as.numeric(cvg.n5[[chr]])
  cvg.n3 <- as.numeric(cvg.n3[[chr]])
  
  pos <- data.frame(seqnames = chr, 
                    start = 1:len, 
                    end = 1:len, 
                    strand = "+", 
                    cvg.pg = cvg.pg, 
                    cvg.p5 = cvg.p5, 
                    cvg.p3 = cvg.p3)
  neg <- data.frame(seqnames = chr, 
                    start = 1:len, 
                    end = 1:len, 
                    strand = "-", 
                    cvg.ng = cvg.ng, 
                    cvg.n5 = cvg.n5, 
                    cvg.n3 = cvg.n3)
  
  return(list(pos = pos, neg = neg))
}

localMax <- function(x, se, w, bg, min.se) {
  len <- length(x)
  idx <- which(x >= bg & se >= min.se)
  
  n <- 0
  lm <- 0
  
  for (i in idx) {
    if (i > w & i < len - w) {
      max <- which.max(x[(i - w):(i + w)]) + (i - w - 1)
    } else if (i <= w) {
      max <- which.max(x[1:(i + w)])
    } else {
      max <- which.max(x[(i - w):len]) + (i - w - 1)
    }
    
    if (max == i) {
      n <- n + 1
      lm[n] <- max
    }
  }
  return(lm)
}

flattopRatio <- function(idx, cvg, right = TRUE) {
  w1 <- sapply(idx, function(x) (x - 1):(x + 1))
  
  if (right) {
    w2 <- sapply(idx, function(x) (x - 1):(x + 15))
  } else {
    w2 <- sapply(idx, function(x) (x - 15):(x + 1))
  }
  
  s1 <- 0
  for (i in 1:ncol(w1)) {
    s1[i] <- sum(cvg[w1[, i]])
  }
  
  s2 <- 0
  for (i in 1:ncol(w2)) {
    s2[i] <- sum(cvg[w2[, i]])
  }
  ratio <- s1/s2
  return(ratio)
}

matchEnds <- function(lm, dat, right = TRUE) {
  if (right) {
    ## if (x + 49 > chrlen) mate.win <- sapply(lm, function(x) (x + 14):chrlen) else
    mate.win <- sapply(lm, function(x) (x + 14):(x + 49))
  } else {
    ## if (x - 49 < 1) mate.win <- sapply(lm, function(x) 1:(x - 14)) else
    mate.win <- sapply(lm, function(x) (x - 49):(x - 14))
  }
  
  mate <- 0
  for (i in 1:ncol(mate.win)) {
    mate[i] <- mate.win[1, i] + which.max(dat[mate.win[, i]]) - 1
  }
  return(mate)
}

getPeaks <- function(dat, chr, w = 15, bg = 20, min.se = 0.5, min.ft = 0.8, Edges = NULL) {
  cvg.p5 <- dat[["pos"]]$cvg.p5
  cvg.p3 <- dat[["pos"]]$cvg.p3
  cvg.n5 <- dat[["neg"]]$cvg.n5
  cvg.n3 <- dat[["neg"]]$cvg.n3
  
  cvg.pg <- dat[["pos"]]$cvg.pg
  cvg.ng <- dat[["neg"]]$cvg.ng
  
  se.p5 <- cvg.p5/cvg.pg
  se.p3 <- cvg.p3/cvg.pg
  se.n5 <- cvg.n5/cvg.ng
  se.n3 <- cvg.n3/cvg.ng
  
  lm.p5 <- localMax(x = cvg.p5, se = se.p5, w = w, bg = bg, min.se = min.se)
  lm.p3 <- localMax(x = cvg.p3, se = se.p3, w = w, bg = bg, min.se = min.se)
  lm.n5 <- localMax(x = cvg.n5, se = se.n5, w = w, bg = bg, min.se = min.se)
  lm.n3 <- localMax(x = cvg.n3, se = se.n3, w = w, bg = bg, min.se = min.se)
  
  ft.p5 <- flattopRatio(idx = lm.p5, cvg = cvg.p5, right = TRUE)
  ft.p3 <- flattopRatio(idx = lm.p3, cvg = cvg.p3, right = FALSE)
  ft.n5 <- flattopRatio(idx = lm.n5, cvg = cvg.n5, right = FALSE)
  ft.n3 <- flattopRatio(idx = lm.n3, cvg = cvg.n3, right = TRUE)
  
  p5 <- lm.p5[which(ft.p5 >= min.ft)]
  p3 <- lm.p3[which(ft.p3 >= min.ft)]
  n5 <- lm.n5[which(ft.n5 >= min.ft)]
  n3 <- lm.n3[which(ft.n3 >= min.ft)]
  
  p5m <- matchEnds(p5, cvg.p3, right = TRUE)
  p3m <- matchEnds(p3, cvg.p5, right = FALSE)
  n5m <- matchEnds(n5, cvg.n3, right = FALSE)
  n3m <- matchEnds(n3, cvg.n5, right = TRUE)
  
  gr.p5 <- GRanges(seqnames = chr, IRanges(start = p5, end = p5m), strand = "+")
  gr.p3 <- GRanges(seqnames = chr, IRanges(start = p3m, end = p3), strand = "+")
  gr.n5 <- GRanges(seqnames = chr, IRanges(start = n5m, end = n5), strand = "-")
  gr.n3 <- GRanges(seqnames = chr, IRanges(start = n3, end = n3m), strand = "-")
  
  if (Edges == 1) {
    gr <- c(gr.p5, gr.n5)
  } else if (Edges == 2) {
    gr <- c(gr.p3, gr.n3)
  } else if (Edges == 3) {
    gr <- unique(c(gr.p5, gr.p3, gr.n5, gr.n3))
  } else if (Edges == 4) {
    ## merge overlap footprints, in case there are.
    ov.p <- findOverlaps(gr.p5, gr.p3, minoverlap = 15L)
    ov.n <- findOverlaps(gr.n5, gr.n3, minoverlap = 15L)
    gr <- unique(c(gr.p5[queryHits(ov.p)], gr.n5[queryHits(ov.n)]))
  } else {
    stop("Error: Parameter 'Edges = ", Edges, "' is not recognized.")
  }
}

# getCoveragesBam(file, selseqlen){
#   seqname <- names(selseqlen)
#   seqlen <- unname(selseqlen)
#   param <- ScanBamParam(what = c("pos", "qwidth", "strand"),
#                         which = GRanges(seqname, IRanges(1,seqlen)),
#                         flag = scanBamFlag(isUnmappedQuery = FALSE))
#   x <- scanBam(file,param = param)[[1]]
#   coverage(IRanges(x[["pos"]], width = x[["qwidth"]]))
# }
