
load_map_file <- function(filename) {
  reads <- read.table(filename, sep="\t", 
                      col.names=c("read.name", "strand", "ref.name", 
                                  "pos", "sequence", "quality", 
                                  "duplicity", "mismatch"),
                      stringsAsFactors=FALSE,
                      colClasses=c("character","character","character",
                                   "integer", "character","character",
                                   "integer","character"))
  reads$reads <- as.integer(vapply(strsplit(reads$read.name, '-', fixed=T), 
                                   function(x) x[2], ""))
  reads$pos <- reads$pos + 1   # bowtie positions are zero-indexed
  return(reads)
}

load_map_files <- function(samples, path) {
  
}
