
library(ShortRead)
library(Biostrings)

filename <- "~/seqdata/strain_comparison/input/kzz_lane1_1.fastq.gz"

sampler <- FastqSampler(filename, 20000)
reads <- yield(sampler)

trim1 <- narrow(reads, start=10, end=34)
trim2 <- trimLRPatterns(Rpattern=DNAString("TAACAG"), subject=trim1, max.Rmismatch=c(-1,-1,1,1,1,1))
unmatched <- width(trim2) > 21
n_unmatched <- sum(as.integer(unmatched))
trim3 <- trim2[!unmatched]
barcode_strings <- sread(narrow(trim3), start=1, end=6)

match_barcode <- function(dnas, barcodes, mismatches=0) {
  ncodes <- length(barcodes)
  ndnas <- length(dnas)
  candidates <- subseq(dnas, start=1, end=nchar(barcodes[1]))
  hits <- matrix(FALSE, nrow=ndnas, ncol=ncodes, dimnames=list(NULL, barcodes))
  for (code in barcodes) {
    hits[ ,code] <- isMatchingStartingAt(code, candidates, starting.at=1, 
                                         max.mismatch=mismatches)
  }
  if (any(rowSums(hits) > 1)) {
    warning("Mismatch level creates ambiguous matches.")
  }
  
  return(hits)
}

barcodes <- c("ACGAAT", "TGCTGG", "GAACTA")
hits <-  match_barcode(sread(trim3), barcodes, mismatches=2)

