
require(stringr)
require(dplyr)

tigr4_19f <- read.delim("~/seqdata/homology/tigr4_19f.txt", stringsAsFactors=F)
names(tigr4_19f) <- c("tigr4", "SPTaiwan_19F")

tigr4_d39 <- read.delim("~/seqdata/homology/tigr4_d39.txt", 
                        stringsAsFactors=F)[,c(2,3)]
names(tigr4_d39) <- c("tigr4", "d39")
tigr4_d39$d39 <- str_trim(tigr4_d39$d39)
tigr4_d39[!str_detect(tigr4_d39$tigr4,"^SP_"),"tigr4"] <- NA
tigr4_d39[!str_detect(tigr4_d39$d39,"^SPD"),"d39"] <- NA

strep_homologs <- full_join(tigr4_19f, tigr4_d39)
set_nas <- function(x) ifelse(is.na(x) | !str_detect(x, "SP"), NA, x)
strep_homologs <- as.data.frame(lapply(strep_homologs, set_nas), 
                                stringsAsFactors=F)
strep_homologs <- strep_homologs[sapply(1:nrow(strep_homologs), 
                                 function(i) !all(is.na(strep_homologs[i,]))),]

make_mapping <- function(homologs) {
  strains <- names(homologs)
  maxlen <- sapply(homologs, function(x) max(sapply(x[!is.na(x)], nchar)))
  pads <- sapply(maxlen, function(x) paste(rep(" ", times=x), collapse=""))
  for (strain in strains) {
    homologs[[strain]][is.na(homologs[[strain]])] <- pads[strain]
  }
  joined <- do.call(paste, as.list(homologs))
  get_key_frame <- function(strain) {
    frame <- data.frame(gene=homologs[[strain]], name=joined, 
                        stringsAsFactors=F)
    frame <- frame[!is.na(frame$gene),]
    return(frame)
  }
  maps <- do.call(rbind, lapply(strains, get_key_frame))
  mapping <- maps$name
  names(mapping) <- maps$gene
  return(mapping)
}

add_homologs <- function(df, homologs) {
  mapping <- make_mapping(homologs)
  df$homolog <- mapping[df$locus]
  return(df)
}

