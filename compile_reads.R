
load.map.file <- function(filename) {
  reads <- read.table(filename, sep="\t", 
                      col.names=c("read.name", "strand", "ref.name", "pos", 
                                  "sequence", "quality", "duplicity", "mismatch"),
                      stringsAsFactors=FALSE,
                      colClasses=c("character","character","character","integer",
                                   "character","character","integer","character"))
  reads$reads <- as.integer(vapply(strsplit(reads$read.name, '-', fixed=T), function(x) x[2], ""))
  reads$pos <- reads$pos + 1   # bowtie positions are zero-indexed
  reads
}

map.reads <- function(reads,genes) {
  # this method assumes that no genes overlap, which should be true
  find.locus <- function(pos) {
    locus <- genes$locus_tag[genes$start < pos & pos < genes$stop]
    if (length(locus) < 1)
      ""
    else
      locus[length(locus)]  # choose the last gene; it is likely in the tail of the other
  }
  loci <- vapply(reads$pos, find.locus, "")
  map <- data.frame(gene=loci, reads=reads$reads, strand=reads$strand, position=reads$pos, stringsAsFactors=FALSE)
  map[vapply(map$gene,nchar,0) > 0,]
}

consolidate.reads <- function(map) ddply(map,.(gene),summarize, reads = sum(reads))

#reads <- load.reads.file("test/L1_2394_WT_FRUC_T1.map")
#map <- map.reads(reads,genes)


