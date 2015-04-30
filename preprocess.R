
JULIA_SPLITTER_CMD <- "julia ~/Dropbox/bc/fastx_tools/tnseq_barcode_splitter.jl"

BOWTIE_CMD <- "~/bowtie/bowtie"
BOWTIE_OPTS <- "-f -m 1 -n 1 --best -y -p 2"
INDEX_PATH <- "~/seqdata/index/"

create_tnseq_experiment <- function(path) {
  # ensure path ends with an "/"
  if (!str_detect(path, "/$")) {
    path <- paste0(path, "/")
  }
  
  # create a directory inside the path, unless it already exists
  create_dir <- function(dirname) {
    fullpath <- paste0(path, dirname, "/")
    if (!file.exists(fullpath)) {
      system(paste("mkdir", fullpath))
    }
    return(path.expand(fullpath))
  }
  
  tnseq <- list(
    path = path,
    input_path = create_dir("input"),
    split_path = create_dir("split"),
    map_path = create_dir("map"),
    log_path = create_dir("log")
  )
  
  # check that barcode file exists in main directory
  barcode_file <- paste0(path, "barcodes.txt")
  if (!file.exists(barcode_file)) {
    stop("Cannot find barcode file: ", barcode_file)
  }
  tnseq$barcode_file <- barcode_file
  
  # check that sample sheet exists in the main directory
  sample_file <- paste0(path, "samples.csv")
  if (!file.exists(sample_file)) {
    stop("Cannot find sample file: ", sample_file)
  }
  tnseq$sample_file <- sample_file
  tnseq$samples <- read.csv(sample_file, stringsAsFactors=F)
  # allow case insensitive names in samples file
  names(tnseq$samples) <- tolower(names(tnseq$samples))
  
  # read all input files and create default filelog
  lane_files <- system(paste("ls -1", tnseq$input_path), intern=T)
  filelog <- data.frame(lane=lane_files, 
                        file=paste0(tnseq$input_path, lane_files),
                        stringsAsFactors=F)
  tnseq$filelog <- list(input=filelog)
  
  return(tnseq)
}


julia_barcode_splitter <- function(tnseq) {
  filelog <- tnseq$filelog$input
  outlog <- NULL
  barcodes <- read.table(tnseq$barcode_file, 
                         col.names=c("name","code"), 
                         as.is=TRUE)$name
  for (i in 1:nrow(filelog)) {
    # julia cannot handle tilde expansion
    logfile <- paste0(tnseq$log_path, filelog$lane[i], "_split.log")
    
    cmd <- paste(JULIA_SPLITTER_CMD,
                 "-m 1 -p 1", 
                 "--prefix", path.expand(filelog$file[i]),
                 "--barcode", path.expand(tnseq$barcode_file), 
                 "-i", path.expand(filelog$file[i]),
                 ">", logfile)
    system(cmd)
    
    newlog <- filelog[rep.int(i,length(barcodes)),]
    newlog$file <- paste0(filelog$file[i], "_", barcodes, ".fastq")
    newlog$barcode <- barcodes
    
    logs <- read.table(logfile, stringsAsFactors=F, header=T)
    reads <- logs$reads
    names(reads) <- logs$pool
    newlog$reads <- reads[newlog$barcode]
    
    outlog <- rbind(outlog, newlog)
  }
  
  tnseq$filelog$debarcoded <- outlog
  return(tnseq)
}

count_reads <- function(filelog) {
  countf <- function(file) {
    output <- system(paste("wc -l", file), intern=T)
    as.integer(str_split(str_trim(output[1]), " ")[[1]][1]) / 4
  }
  filelog$reads <- sapply(filelog$file, countf)
  return(filelog)
}

cleanup_files <- function(filelog) {
  sapply(filelog$file, file.remove)
}

final_filter_collapse <- function(tnseq) {
  filelog <- tnseq$filelog$debarcoded
  quallog <- filelog
  quallog$file <- paste0(filelog$file, ".quality")
  for (i in 1:nrow(filelog)) {
    if (filelog$reads[i] > 0) {
      cmd <- paste("(",
                   "fastq_quality_filter -Q 33 -q 8 -p 100", "|",
                   "fastx_collapser -Q 33",
                   ")",
                   "<", filelog$file[i], 
                   ">", quallog$file[i])
    } else {
      cmd <- paste("cp", filelog$file[i], quallog$file[i])
    }
    #cat(cmd)
    system(cmd)
  }
  cleanup_files(filelog)
  quallog$reads <- filelog$reads
  
  tnseq$filelog$filtered_collapsed <- quallog
  return(tnseq)
}

move_split_files <- function(tnseq) {
  filelog <- tnseq$filelog$filtered_collapsed
  splitpath <- tnseq$split_path
  
  newfile <- character(nrow(filelog))
  for (i in 1:nrow(filelog)) {
    newfile[i] <- paste0(splitpath,
                         paste(filelog$lane[i],  
                               filelog$barcode[i], sep="_"),
                         ".fastq")
    system(paste("mv", filelog$file[i], newfile[i]))
  }
  newlog <- filelog
  newlog$file <- newfile
  
  tnseq$filelog$split <- newlog
  return(tnseq)
}


preprocess_reads <- function(tnseq) {
  if (nrow(tnseq$filelog$input) > 1) {
    # iterate over one input file at a time
    # first, expand into multiple tnseq objects, each with a single input
    n <- nrow(tnseq$filelog$input)
    tnseqs <- lapply(1:n, function(x) tnseq)
    for (i in 1:n) {
      tnseqs[[i]]$filelog$input <- tnseqs[[i]]$filelog$input[i,]
    }
    
    # run preprocess_reads on each input file
    tnseqs <- lapply(tnseqs, preprocess_reads)
    
    # combine the results back into a single tnseq object
    tnseq <- tnseqs[[1]]
    for (name in names(tnseq$filelog)) {
      dfs <- lapply(tnseqs, function(x) x$filelog[[name]])
      tnseq$filelog[[name]] <- do.call(rbind, dfs)
    }
  } else {
    # preprocess a single input file
    print(paste("Processing file", tnseq$filelog$input$file[1]))
    ptm <- proc.time()
    tnseq <- julia_barcode_splitter(tnseq)
    tnseq <- final_filter_collapse(tnseq)
    tnseq <- move_split_files(tnseq)
    print(proc.time() - ptm)
  }

  return(tnseq)
}


# ================ sample mapping ================

parse_bowtie_log <- function(logfile) {
  text <- paste(readLines(logfile), collapse="")
  reads <- as.integer(str_match(text, "reads processed: (\\d+)")[1,2])
  aligned <- as.integer(str_match(text, "alignment: (\\d+)")[1,2])
  return(c(reads=reads, aligned=aligned))
}

map_reads <- function(tnseq) {
  n_samples <- nrow(tnseq$samples)
  tnseq$samples$mapfile <- character(n_samples)
  tnseq$samples$reads <- integer(n_samples)
  tnseq$samples$aligned <- integer(n_samples)
  for (i in 1:n_samples) {
    lane <- tnseq$samples$lane[i]
    code <- tnseq$samples$barcode[i]
    split_idx <- which(tnseq$filelog$split$lane == lane & 
                         tnseq$filelog$split$barcode == code)
    if (length(split_idx) > 1) {
      stop("Multiple samples for lane ", lane, " barcode ", code)
    }
    splitfile <- tnseq$filelog$split$file[split_idx]
    mapfile <- paste0(tnseq$map_path, lane, "_", code, ".map")
    logfile <- paste0(tnseq$log_path, lane, "_", code, "_bowtie.log")
    indexfile <- paste0(INDEX_PATH, tnseq$samples$genome[i])
    system(paste(BOWTIE_CMD, BOWTIE_OPTS, indexfile, splitfile, mapfile,
                 "2>", logfile))
    
    tnseq$samples$mapfile[i] <- mapfile
    results <- parse_bowtie_log(logfile)
    tnseq$samples$reads[i] <- results["reads"]
    tnseq$samples$aligned[i] <- results["aligned"]
  }
  return(tnseq)
}


#path <- "~/seqdata/strain_comparison"
#path <- "~/Dropbox/bc/tnseqr/test"
#path <- "~/seqdata/strain_comparison_nextseq_3_31"

#tnseq <- create_tnseq_experiment(path)

#ptm <- proc.time()
#tnseq <- preprocess_reads(tnseq)
#elapsed <- proc.time() - ptm

#tnseq <- map_reads(tnseq)

