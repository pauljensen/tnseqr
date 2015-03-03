
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
  }
  file.remove(paste0(filelog$file[1], ".unmatched"))
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

move_split_files <- function(filelog, path) {
  splitpath <- paste0(path, "/split")
  system(paste("mkdir", splitpath))
  newfile <- character(nrow(filelog))
  for (i in 1:nrow(filelog)) {
    newfile[i] <- paste0(splitpath, "/",
                         paste(filelog$lane[i], 
                               filelog$overhang[i], 
                               filelog$barcode[i], sep="_"),
                         ".fastq")
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
                       file=paste(path, "input", lane_files, sep="/"))
  return(filelog)
}

filelog1 <- preprocess("test")
filelog2 <- trim_front_back(filelog1)
filelog3 <- split_by_overhang_length(filelog2)
filelog4 <- trim_transposon_sequence(filelog3)
filelog5 <- split_by_barcode(filelog4, "test/barcodes.txt")
filelog6 <- final_trim_filter(filelog5)
filelog7 <- move_split_files(filelog6,"test")
filelog8 <- collapse_by_overhang(filelog7,"test")
