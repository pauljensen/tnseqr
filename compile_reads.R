
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

# ==================== working with MAP files ====================

# load a full map file, output from bowtie
load_map_file <- function(filename) {
  reads <- read.table(filename, sep="\t", 
                      col.names=c("read_name", "strand", "ref_name", 
                                  "pos", "sequence", "quality", 
                                  "duplicity", "mismatch"),
                      stringsAsFactors=FALSE,
                      colClasses=c("character", "factor", "character",
                                   "integer", "character", "character",
                                   "integer","character"))
  reads$reads <- as.integer(str_match(reads$read_name, "\\d+-(\\d+)")[,2])
  reads$pos <- reads$pos + 1   # bowtie positions are zero-indexed
  return(reads)
}

compress_reads <- function(reads) {
  reads %>% 
    select(pos, reads, strand) %>% 
    group_by(pos, strand) %>% 
    summarize(reads=sum(reads))
}

compress_map_file <- function(filename) {
  outfile <- paste0(filename, ".compressed")
  reads <- load_map_file(filename)
  reads$pos <- as.character(reads$pos)
  write.csv(compress_reads(reads), file=outfile, row.names=F, quote=F)
}

compress_map_files <- function(tnseq) {
  mclapply(tnseq$samples$mapfile, compress_map_file)
  return(tnseq)
}

load_compressed_map_file <- function(filename) {
  reads <- read.csv(filename, header=T, 
                    col.names=c("pos", "strand", "reads"),
                    colClasses=c("character", "factor", "integer"))
  reads$pos <- as.integer(reads$pos)
  return(reads)
}

get_read_pairs <- function(mapfile1, mapfile2) {
  compressed1 <- paste0(mapfile1, ".compressed")
  compressed2 <- paste0(mapfile2, ".compressed")
  if (file.exists(compressed1)) {
    reads1 <- load_compressed_map_file(compressed1)
  } else {
    reads1 <- load_map_file(mapfile1) %>% compress_reads()
  }
  if (file.exists(compressed2)) {
    reads2 <- load_compressed_map_file(compressed2)
  } else {
    reads2 <- load_map_file(mapfile2) %>% compress_reads()
  }
  
  both <- full_join(reads1, reads2, by=c("pos", "strand"))
  both$reads1 <- both$reads.x
  both$reads.x <- NULL
  both$reads2 <- both$reads.y
  both$reads.y <- NULL

  both$reads1[is.na(both$reads1)] <- 0
  both$reads2[is.na(both$reads2)] <- 0
  
  return(both)
}

# =================== working with individual insertions =====================

load_insertions <- function(tnseq) {
  n_paired <- nrow(tnseq$paired_samples)
  cols_to_keep <- c("strain", "library", "condition", "genome", "expansion")
  load_aux <- function(i) {
    cat(i, "\n")
    inserts <- get_read_pairs(tnseq$paired_samples$t1_mapfile[i],
                              tnseq$paired_samples$t2_mapfile[i])
    df <- do.call(data.frame, tnseq$paired_samples[i,cols_to_keep])
    cbind(df, inserts)
  }
  tnseq$insertions <- do.call(rbind, mclapply(1:nrow(tnseq$paired_samples), 
                                              load_aux))
  return(tnseq)
}

filter_insertions <- function(tnseq, drop_t1_zeros=T, drop_t2_zeros=F, min_total_reads=35) {
  to_drop <- !logical(nrow(tnseq$insertions))
  if (drop_t1_zeros) {
    to_drop <- to_drop & (tnseq$insertions$reads1 == 0)
  }
  if (drop_t2_zeros) {
    to_drop <- to_drop & (tnseq$insertions$read2 == 0)
  }
  if (min_total_reads > 0) {
    to_drop <- to_drop & (tnseq$insertions$reads1 + tnseq$insertions$reads2 < 
                          min_total_reads)
  }
  tnseq$insertions <- tnseq$insertions[!to_drop,]
  return(tnseq)
}

map_insertions_to_genes <- function(tnseq, ignore_start=0.0, ignore_end=0.1) {
  genomes <- unique(tnseq$paired_samples$genome)
  genome_files <- paste0("~/seqdata/genomes/", genomes, ".gbk")
  genes <- lapply(genome_files, load_genbank_file)
  names(genes) <- genomes
  
  genome2irange <- function(g) {
    st <- g$start
    en <- g$stop
    st2 <- st
    st[en < st2] <- en[en < st2]
    en[en < st2] <- st2[en < st2]
    width <- en - st
    en <- en - ceiling(ignore_end * width)
    st <- st + floor(ignore_start * width)
    return(IRanges(start=st, end=en, names=g$locus))
  }
  
  iranges <- lapply(genes, genome2irange)
  
  map_aux <- function(df) {
    locs <- IRanges(start=df$pos, width=1)
    ov <- findOverlaps(locs, iranges[[df$genome[1]]])
    n <- nrow(df)
    df$locus <- character(n)
    df$locus[] <- NA
    df$start_pos <- integer(n)
    df$start_pos[] <- NA
    df$end_pos <- df$start_pos
    df$frac <- as.numeric(df$start_pos)
    query <- queryHits(ov)
    subject <- subjectHits(ov)
    df$locus[query] <- as.character(genes[[df$genome[1]]]$locus[subject])
    df$start_pos[query] <- genes[[df$genome[1]]]$start[subject]
    df$end_pos[query] <- genes[[df$genome[1]]]$stop[subject]
    df$frac <- (df$pos - df$start_pos) / (df$end_pos - df$start_pos)
    return(df)
  }
  
  ptm <- proc.time()
  tnseq$insertions %<>%
    group_by(genome) %>%
    do(map_aux(.)) %>%
    ungroup()
  print(proc.time() - ptm)
  
  return(tnseq)
}

calculate_fitness <- function(tnseq, group_libraries=T) {
  tnseq$fitness <- tnseq$insertions %>%
    group_by_("strain", "library", "condition") %>%
    mutate(f1=reads1/sum(reads1),
           f2=reads2/sum(reads2),
           W=ifelse(reads2==0, 0,log(mean(expansion) * f2/f1) / 
                                   log(mean(expansion) * (1-f2)/(1-f1)))) %>%
    ungroup()
  tnseq$fitness %<>% 
    group_by_("strain", "library", "condition", "locus") %>%
    summarize(
      insertions=n(),
      fitness=mean(W),
      fitness_sd=sd(W)
    ) %>% ungroup()

  return(tnseq)
}

#reads <- load.reads.file("test/L1_2394_WT_FRUC_T1.map")
#map <- map.reads(reads,genes)

#r1 <- "~/seqdata/strain_comparison/map/kzz_lane1.fastq_BCINDX1a.map"
#r2 <- "~/seqdata/strain_comparison/map/kzz_lane1.fastq_BCINDX2c.map"
#reads <- get_read_pairs(r1, r2)
#genes <- load_genbank_file("~/seqdata/genomes/tigr4.gbk")
#fit <- calculate_gene_fitness(reads, genes, expansion=1000)

