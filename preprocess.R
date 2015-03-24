
library(stringr)
library(dplyr)

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

fastx_trimmer <- function(Q=NA, f=NA, t=NA, i=NA, l=NA) {
  fxopt <- function(x) {
    quoted <- substitute(x)
    if (!is.na(x))
      paste0("-", deparse(quoted), " ", x)
    else
      ""
  }
  paste("fastx_trimmer", fxopt(Q), fxopt(f), fxopt(t), fxopt(i), fxopt(l))
}

trim_front_back <- function(filelog) {
  trimfiles <- paste0(filelog$file, ".trimmed")
  for (i in 1:length(filelog$file)) {
    cmd <- paste(fastx_trimmer(Q=33, f=10, i=filelog$file[i]), "|",
                 fastx_trimmer(Q=33, t=17), ">", trimfiles[i])
    system(cmd)
  }
  filelog$file <- trimfiles
  return(count_reads(filelog))
}

split_by_overhang_length <- function(filelog) {
  lengths <- as.character(14:17)
  
  barcodes <- c("TAACAG", "TAACA", "TAAC", "TAA")
  names(barcodes) <- lengths
  
  outlog <- NULL
  for (i in 1:nrow(filelog)) {
    for (len in lengths) {
      system(paste0('echo ', len, "\t", barcodes[len], " > barcode.txt"))
      cmd <- paste("fastx_barcode_splitter.pl",
                   "--eol",
                   "--mismatches 1",
                   "--prefix", paste0(filelog$file[i], "."),
                   "--bcfile barcode.txt",
                   "<", filelog$file[i])
      system(cmd)
    }
    newlog <- filelog[rep.int(i,length(lengths)),]
    newlog$file <- paste0(filelog$file[i], ".", lengths)
    newlog$overhang <- lengths
    outlog <- rbind(outlog, newlog)
    file.remove(paste0(filelog$file[i], ".unmatched"))
  }
  file.remove("barcode.txt")
  cleanup_files(filelog)
  return(count_reads(outlog))
}

trim_transposon_sequence <- function(filelog) {
  trimmed <- paste0(filelog$file, ".trimmed")
  for (i in 1:nrow(filelog)) {
    len <- as.integer(filelog$overhang[i])
    cmd <- paste(fastx_trimmer(Q=33, l=as.integer(len)+7),
                 "<", filelog$file[i],
                 ">", trimmed[i])
    system(cmd)
  }
  newlog <- filelog
  newlog$file <- trimmed
  cleanup_files(filelog)
  return(newlog)
}

load_barcodes <- function(bcfile) {
  read.table(bcfile, col.names=c("name","code"), as.is=TRUE)$name
}

split_by_barcode <- function(filelog, bcfile) {
  outlog <- NULL
  barcodes <- load_barcodes(bcfile)
  for (i in 1:nrow(filelog)) {
    cmd <- paste("fastx_barcode_splitter.pl --bol --mismatches 1 --partial 1", 
                 "--prefix", paste0(filelog$file[i], "."),
                 "--bcfile", bcfile, 
                 "<", filelog$file[i])
    system(cmd)
    newlog <- filelog[rep.int(i,length(barcodes)),]
    newlog$file <- paste0(filelog$file[i], ".", barcodes)
    newlog$barcode <- barcodes
    outlog <- rbind(outlog, newlog)
    file.remove(paste0(filelog$file[i], ".unmatched"))
  }
  cleanup_files(filelog)
  return(count_reads(outlog))
}

final_trim_filter <- function(filelog) {
  trimlog <- filelog
  trimlog$file <- paste0(filelog$file, ".trimmed")
  for (i in 1:nrow(filelog)) {
    if (filelog$reads[i] > 0) {
      cmd <- paste(fastx_trimmer(Q=33, f=7), "<", filelog$file[i], 
                   ">", trimlog$file[i])
    } else {
      cmd <- paste("cp", filelog$file[i], trimlog$file[i])
    }
    system(cmd)
  }
  trimlog <- count_reads(trimlog)
  cleanup_files(filelog)
  
  quality <- paste0(trimlog$file, ".quality")
  for (i in 1:nrow(trimlog)) {
    if (filelog$reads[i] > 0) {
      cmd <- paste("fastq_quality_filter -Q 33 -q 8 -p 100 -o", quality[i],
                   "<", trimlog$file[i])
    } else {
      cmd <- paste("cp", trimlog$file[i], quality[i])
    }
    system(cmd)
  }
  cleanup_files(trimlog)
  quallog <- trimlog
  quallog$file <- quality
  filelog <- count_reads(quallog)
  filelog$quality_reads <- filelog$reads
  filelog$reads <- quallog$reads
  return(filelog)
}

final_trim_filter_collapse <- function(filelog) {
  quallog <- filelog
  quallog$file <- paste0(filelog$file, ".quality")
  for (i in 1:nrow(filelog)) {
    if (filelog$reads[i] > 0) {
      cmd <- paste("(",
                   fastx_trimmer(Q=33, f=7), "|",
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
  quallog$quality_reads <- quallog$reads
  quallog$reads <- filelog$reads
  return(quallog)
}

final_filter_collapse <- function(filelog) {
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
  quallog$quality_reads <- quallog$reads
  quallog$reads <- filelog$reads
  return(quallog)
}

move_split_files <- function(filelog, path, overhang=T) {
  splitpath <- paste0(path, "/split")
  system(paste("mkdir", splitpath))
  newfile <- character(nrow(filelog))
  for (i in 1:nrow(filelog)) {
    if (overhang) {
      newfile[i] <- paste0(splitpath, "/",
                           paste(filelog$lane[i], 
                                 filelog$overhang[i], 
                                 filelog$barcode[i], sep="_"),
                           ".fastq")
    } else {
      newfile[i] <- paste0(splitpath, "/",
                           paste(filelog$lane[i],  
                                 filelog$barcode[i], sep="_"),
                           ".fastq")
    }
    system(paste("mv", filelog$file[i], newfile[i]))
  }
  newlog <- filelog
  newlog$file <- newfile
  return(newlog)
}

collapse_by_overhang <- function(filelog, path) {
  newlog <- filelog %>% 
    group_by(lane, barcode) %>%
    summarize(reads=sum(reads), 
              quality_reads=sum(quality_reads),
              catargs=paste(file, collapse=" "))
  
  newfiles <- paste0(path, "/split/", 
                     newlog$lane, "_", newlog$barcode, ".fastq")
  for (i in 1:nrow(newlog)) {
    system(paste("cat", newlog$catargs[i], ">", newfiles[i]))
  }
  newlog$catargs <- NULL
  newlog$file <- newfiles
  cleanup_files(filelog)
  return(as.data.frame(newlog))
}

preprocess <- function(path) {
  lane_files <- system(paste0("ls -1 ", path, "/input"), intern=T)
  filelog <- data.frame(lane=lane_files, 
                       file=paste(path, "input", lane_files, sep="/"),
                       stringsAsFactors=F)
  return(filelog)
}

JULIA_SPLITTER_CMD <- "julia ~/Dropbox/bc/fastx_tools/tnseq_barcode_splitter.jl"

julia_barcode_splitter <- function(filelog, bcfile) {
  outlog <- NULL
  barcodes <- load_barcodes(bcfile)
  for (i in 1:nrow(filelog)) {
    # julia cannot handle tilde expansion
    cmd <- paste(JULIA_SPLITTER_CMD,
                 "-m 1 -p 1", 
                 "--prefix", path.expand(filelog$file[i]),
                 "--barcode", path.expand(bcfile), 
                 "-i", path.expand(filelog$file[i]))
    system(cmd)
    
    newlog <- filelog[rep.int(i,length(barcodes)),]
    newlog$file <- paste0(filelog$file[i], "_", barcodes, ".fastq")
    newlog$barcode <- barcodes
    
    #logs <- read.table("out.log",stringsAsFactors=F)
    #reads <- logs$reads
    #names(reads) <- logs$pool
    #cat(reads)
    #newlog$reads <- reads[newlog$barcode]
    #file.remove("out.log")
    
    outlog <- rbind(outlog, newlog)
    #file.remove(paste0(filelog$file[i], ".unmatched"))
  }
  #cleanup_files(filelog)
  return(count_reads(outlog))
}

process_tn_seq <- function(path) {
  bcfile <- paste0(path, "/barcodes.txt")
  
  filelog <- vector("list", length=8)
  filelog[[1]] <- preprocess(path)
  filelog[[2]] <- trim_front_back(filelog[[1]])
  filelog[[3]] <- split_by_overhang_length(filelog[[2]])
  filelog[[4]] <- trim_transposon_sequence(filelog[[3]])
  filelog[[5]] <- split_by_barcode(filelog[[4]], bcfile)
  filelog[[6]] <- final_trim_filter_collapse(filelog[[5]])
  filelog[[7]] <- move_split_files(filelog[[6]], path)
  filelog[[8]] <- collapse_by_overhang(filelog[[7]], path)
  
  return(filelog)
}

process_tn_seq_fast <- function(path) {
  bcfile <- paste0(path, "/barcodes.txt")
  filelog <- vector("list", length=1)
  filelog[[1]] <- preprocess(path)
  filelog[[2]] <- julia_barcode_splitter(filelog[[1]], bcfile)
  filelog[[3]] <- final_filter_collapse(filelog[[2]])
  filelog[[4]] <- move_split_files(filelog[[3]], path, overhang=F)
  
  return(filelog)
}

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
#samples <- map_reads(load_sample_sheet(path), path)

#filelog <- process_tn_seq_fast("~/Dropbox/bc/tnseqr/test")

ptm <- proc.time()
filelog <- process_tn_seq(path)
elapsed <- proc.time() - ptm

#filelog1 <- preprocess("~/seqdata/strain_comparison")
#filelog2 <- trim_front_back(filelog1)
#filelog3 <- split_by_overhang_length(filelog2)
#filelog4 <- trim_transposon_sequence(filelog3)
#filelog5 <- split_by_barcode(filelog4, "~/seqdata/strain_comparison/barcodes.txt")
#filelog6 <- final_trim_filter_collapse(filelog5)
#filelog7 <- move_split_files(filelog6,"~/seqdata/strain_comparison")
#filelog8 <- collapse_by_overhang(filelog7,"~/seqdata/strain_comparison/test")
