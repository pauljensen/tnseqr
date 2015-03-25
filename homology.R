
tigr4_19f <- read.delim("~/seqdata/homology/tigr4_19f.txt", stringsAsFactors=F)
names(tigr4_19f) <- c("tigr4", "SPTaiwan_19F")

tigr4_d39 <- read.delim("~/seqdata/homology/tigr4_d39.txt", stringsAsFactors=F)[,c(2,3)]
names(tigr4_d39) <- c("tigr4", "d39")
tigr4_d39$d39 <- str_trim(tigr4_d39$d39)
tigr4_d39[!str_detect(tigr4_d39$tigr4,"^SP_"),"tigr4"] <- NA
tigr4_d39[!str_detect(tigr4_d39$d39,"^SPD"),"d39"] <- NA

homologs <- full_join(tigr4_19f, tigr4_d39)
max_lengths <- lapply(homologs, function(x) max(sapply(as.character(x), nchar)))
prefixes <- c(tigr4="SP_", SPTaiwan_19F="SPT_", d39="SPD_")

for (strain in names(homologs)) {
  homologs[is.na(homologs[[strain]]),strain] <- "X"
  for (i in 1:nrow(homologs)) {
    if (!str_detect(homologs[i,strain], prefixes[strain])) {
      homologs[i,strain] <- paste(rep(" ", max_lengths[strain]), collapse="")
    }
  }
}

homologs$homolog <- do.call(paste, homologs)

all_genes <- unique(do.call(c, homologs))

add_homologs <- function(fitness, homologs) {
  n <- nrow(fitness)
  fitness$homolog <- fitness$locus
  for (i in 1:n) {
    idx <- which(fitness$locus[i] == homologs[[fitness$genome[i]]])
    if (length(idx) == 1) {
      fitness$homolog[i] <- homologs$homolog[idx]
    }
  }
  return(fitness)
}
