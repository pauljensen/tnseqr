
library(stringr)
library(pryr)

FILES_TO_CLEANUP = character(0)

add_for_cleanup <- function(filename) {
  FILES_TO_CLEANUP <- c(FILES_TO_CLEANUP, filename)
}

fxopt <- function(x) {
  quoted <- substitute(x)
  if (!is.na(x))
    paste0("-", deparse(quoted), " ", x)
  else
    ""
}

fastx_trimmer <- function(Q=NA, f=NA, t=NA, i=NA, l=NA) {
  paste("fastx_trimmer", fxopt(Q), fxopt(f), fxopt(t), fxopt(i), fxopt(l))
}

trim_front_back <- function(filename) {
  trimfile <- paste0(filename, ".trimmed")
  cmd <- paste(fastx_trimmer(Q=33, f=10, i=filename), "|",
               fastx_trimmer(Q=33, t=17), ">", trimfile)
  system(cmd)
  add_for_cleanup(trimfile)
  return(trimfile)
}

parse_fastx_barcode_splitter_output <- function(output) {
  getn <- f(n,line, str_split(output[line],"\t")[[1]][n])
  getcount <- f(line, as.integer(getn(2, line)))
  n_codes <- length(output) - 3
  counts <- numeric()
  counts["unmatched"] <- getcount(n_codes + 2)
  counts["total"] <- getcount(n_codes + 3)
  for (i in (1:n_codes)+1) {
    counts[getn(1,i)] <- getcount(i)
  }
  return(counts)
}

split_by_overhang_length <- function(filename, bcfiles) {
  lengths <- as.character(14:17)
  
  counts <- integer(length(lengths)+1)
  names(counts) <- c(lengths, "total")
  
  outfiles <- paste0(filename, ".", lengths)
  names(lengths) <- lengths
  for (len in lengths) {
    cmd <- paste("fastx_barcode_splitter.pl",
                 "--eol",
                 "--mismatches 1",
                 "--prefix", paste0(filename, "."),
                 "--bcfile", bcfiles[len],
                 "<", filename)
    output <- system(cmd, intern=TRUE)
    counts[len] <- parse_fastx_barcode_splitter_output(output)[len]
    
    # trim off everything after the transposon
    cmd2 <- paste(fastx_trimmer(Q=33, l=as.integer(len)+7),
                  "<", paste0(filename, ".", len),
                  ">", paste0(filename, ".", len, ".trimmed"))
    system(cmd2)
  }
  counts["total"] <- parse_fastx_barcode_splitter_output(output)["total"]
  return(counts)
}





load_barcodes <- function(filename) {
  read.table(filename, col.names=c("name","code"), as.is=TRUE)$name
}

debarcode2 <- function(bcfile) {
  for (len in c("14", "15", "16", "17")) {
    file <- paste0("TNSEQR-", len, ".txt")
    cmd <- paste("fastx_barcode_splitter.pl --bol --mismatches 1 --partial 1", 
                 "--prefix", paste0("'raw/", len, "/'"),
#                 "--prefix", "raw/",
                 "--bcfile", bcfile, 
                 "<", file)
    system(cmd)
  }
}

final_trim_filter <- function(barcodes) {
  pre <- vector("list",4)
  names(pre) <- as.character(14:17)
  post <- vector("list",4)
  names(post) <- as.character(14:17)
  for (len in as.character(14:17)) {
    paths <- paste0("raw/", len, "/", barcodes)
    pre[[len]] <- file.info(paths)
    for (path in paths) {
      cmd <- paste("fastx_trimmer -Q 33 -f 7 |",
                   "fastq_quality_filter -Q 33 -q 8 -p 100 -o", paste0(path,".quality"),
                   "<", path)
      system(cmd)
    post[[len]] <- file.info(paste0(paths, ".quality"))
    }
  }
  return(list(pre=pre, post=post))
}

preprocess_reads <- function(readfiles, path) {
  
}

trimfile <- trim_front_back("test/small-01.txt")
bcfiles <- paste0("bc", as.character(14:17), ".txt")
names(bcfiles) <- as.character(14:17)
counts <- split_by_overhang_length(trimfile, bcfiles)
#debarcode1()
# calculated percent unmapped
#debarcode2("test/bc-01.txt")
#fd <- final_trim_filter(load_barcodes("test/bc-01.txt"))

