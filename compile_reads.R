
library(stringr)
library(dplyr)

pair_samples <- function(tnseq) {
  t1 <- tnseq$samples[tnseq$samples$time == 1,]
  t2 <- tnseq$samples[tnseq$samples$time == 2,]
  n2 <- nrow(t2)
  
  pairs <- t2
  pairs$lane <- NULL
  pairs$barcode <- NULL
  pairs$time <- NULL
  pairs$t2_mapfile <- pairs$mapfile
  pairs$t1_mapfile <- character(n2)
  pairs$mapfile <- NULL
  pairs$t2_reads <- pairs$reads
  pairs$t1_reads <- integer(n2)
  pairs$reads <- NULL
  pairs$t2_aligned <- pairs$aligned
  pairs$t1_aligned <- integer(n2)
  pairs$aligned <- NULL
  
  t1_sl <- paste(t1$strain, t1$library)
  t2_sl <- paste(t2$strain, t2$library)
  t1_slc <- paste(t1_sl, t1$condition)
  t2_slc <- paste(t2_sl, t2$condition)
  
  for (i in 1:n2) {
    idx <- which(t2_slc[i] == t1_slc)
    if (length(idx) == 0) {
      idx <- which(t2_sl[i] == t1_sl)
    }
    if (length(idx) == 0) {
      stop("No matching time 1 for sample: ", t2_slc[i])
    }
    if (length(idx) > 1) {
      stop("Multiple matching time 1 for sample: ", t2_slc[i])
    }
    pairs$t1_mapfile[i] <- t1$mapfile[idx]
    pairs$t1_reads[i] <- t1$reads[idx]
    pairs$t1_aligned[i] <- t1$reads[idx]
  }
  
  tnseq$paired_samples <- pairs
  return(tnseq)
}

load_map_file <- function(filename) {
  reads <- read.table(filename, sep="\t", 
                      col.names=c("read_name", "strand", "ref_name", 
                                  "pos", "sequence", "quality", 
                                  "duplicity", "mismatch"),
                      stringsAsFactors=FALSE,
                      colClasses=c("character", "character", "character",
                                   "integer", "character", "character",
                                   "integer","character"))
  reads$reads <- as.integer(str_match(reads$read_name, "\\d+-(\\d+)")[,2])
  reads$pos <- reads$pos + 1   # bowtie positions are zero-indexed
  return(reads)
}

get_read_pairs <- function(mapfile1, mapfile2) {
  reads1 <- load_map_file(mapfile1)[,c("pos", "reads", "strand")]
  reads2 <- load_map_file(mapfile2)[,c("pos", "reads", "strand")]
  reads1$signed <- reads1$pos * ifelse(reads1$strand == "+", 1, -1)
  reads2$signed <- reads2$pos * ifelse(reads2$strand == "+", 1, -1)
  both <- inner_join(reads1 %>% group_by(pos, strand) %>% 
                       summarize(reads=sum(reads)),
                     reads2 %>% group_by(pos, strand) %>% 
                       summarize(reads=sum(reads)), 
                     by=c("pos", "strand"))
  both$reads1 <- both$reads.x
  both$reads.x <- NULL
  both$reads2 <- both$reads.y
  both$reads.y <- NULL
  return(both)
}

calculate_experiment_fitness <- function(reads, genes, expansion, 
                                         cutoff=0, front_trim=0, end_trim=0.1) {
  f1 <- reads$reads1 / sum(reads$reads1)
  f2 <- reads$reads2 / sum(reads$reads2)
  reads$W <- log(expansion * f2/f1) / log(expansion * (1-f2)/(1-f1))
  loci <- sapply(reads$pos, function(x) which(x >= genes$start & 
                                                x <= genes$stop)[1])
  reads$locus <- genes$locus[loci]
  reads$start <- genes$start[loci]
  reads$stop <- genes$stop[loci]
  reads$frac <- (reads$pos - reads$start) / (reads$stop - reads$start)
  rowfilter <- !is.na(reads$locus) &
                   reads$frac >= front_trim &
                   reads$frac <= 1-end_trim
  reads <- reads[rowfilter,]
  t_test <- function(w) {
    if (length(w) == 1) {
      return(NA)
    } else {
      return(tryCatch(t.test(w, mu=1)$p.value, error=function(e) NA))
    }
  }
  final <- reads %>% group_by(locus) %>% 
              summarise(fitness=mean(W),
                        fitness_sd=sd(W),
                        fitness_p=t_test(W),
                        insertions=n(),
                        reads1=sum(reads1),
                        reads2=sum(reads2))
  return(final)
}

calculate_fitness <- function(tnseq, ...) {
  n_paired <- nrow(tnseq$paired_samples)
  cols_to_keep <- c("strain", "library", "condition", "genome", "expansion")
  
  genomes <- unique(tnseq$paired_samples$genome)
  genome_files <- paste0("~/seqdata/genomes/", genomes, ".gbk")
  genes <- lapply(genome_files, load_genbank_file)
  names(genes) <- genomes
  
  prep_sample <- function(idx) {
    cat("Prepping sample", idx, "\n")
    reads <- get_read_pairs(tnseq$paired_samples$t1_mapfile[idx],
                            tnseq$paired_samples$t2_mapfile[idx])
    expansion <- tnseq$paired_samples$expansion[idx]
    gene <- genes[[tnseq$paired_samples$genome[idx]]]
    fit <- calculate_experiment_fitness(reads, gene, expansion, ...)
    cbind(tnseq$paired_samples[rep.int(idx, nrow(fit)),cols_to_keep], fit)
  }
  tnseq$fitness <- do.call(rbind, lapply(1:n_paired, prep_sample))
  return(tnseq)
}


#reads <- load.reads.file("test/L1_2394_WT_FRUC_T1.map")
#map <- map.reads(reads,genes)

#r1 <- "~/seqdata/strain_comparison/map/kzz_lane1.fastq_BCINDX1a.map"
#r2 <- "~/seqdata/strain_comparison/map/kzz_lane1.fastq_BCINDX2c.map"
#reads <- get_read_pairs(r1, r2)
#genes <- load_genbank_file("~/seqdata/genomes/tigr4.gbk")
#fit <- calculate_gene_fitness(reads, genes, expansion=1000)

