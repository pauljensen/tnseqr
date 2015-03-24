
library(stringr)
library(dplyr)


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
  
  # read all input files and create default filelog
  lane_files <- system(paste("ls -1", tnseq$input_path), intern=T)
  filelog <- data.frame(lane=lane_files, 
                        file=paste0(tnseq$input_path, lane_files),
                        stringsAsFactors=F)
  tnseq$filelog <- list(input=filelog)
  
  return(tnseq)
}


JULIA_SPLITTER_CMD <- "julia ~/Dropbox/bc/fastx_tools/tnseq_barcode_splitter.jl"

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
  tnseq <- julia_barcode_splitter(tnseq)
  tnseq <- final_filter_collapse(tnseq)
  tnseq <- move_split_files(tnseq)
  
  return(tnseq)
}


# ================ sample mapping ================

load_sample_sheet <- function(path) {
  read.csv(paste0(path, "/samples.csv"), stringsAsFactors=F)
}

map_reads <- function(samples, path) {
  map_path <- paste0(path, "/mapped")
  samples$mapfile <- character(nrow(samples))
  system(paste("mkdir", map_path))
  for (i in 1:nrow(samples)) {
    split_file <- paste0(path, "/split/", samples$Lane[i], 
                         "_", samples$Barcode[i], ".fastq")
    map_file <- paste0(map_path, "/", samples$Lane[i], "_",
                       samples$Barcode[i], ".map")
    bowtie_opts <- "-f -m 1 -n 1 --best -y -p 2"
    index_path <- paste0("~/seqdata/index/", samples$Genome[i])
    bowtie_path <- "~/bowtie/bowtie"
    system(paste(bowtie_path, bowtie_opts, index_path,
                 split_file, ">", map_file))
    samples$mapfile[i] <- map_file
  }
  return(samples)
}



path <- "~/seqdata/strain_comparison_test"
#path <- "~/Dropbox/bc/tnseqr/test"

tnseq <- create_tnseq_experiment(path)

#samples <- map_reads(load_sample_sheet(path), path)
#filelog <- process_tn_seq_fast("~/Dropbox/bc/tnseqr/test")

ptm <- proc.time()
tnseq <- preprocess_reads(tnseq)
elapsed <- proc.time() - ptm

