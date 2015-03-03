# genes is from parse_genbank.R
reads1 <- load.map.file("test/0713_L1_T4_WT_T1.map")
map1 <- map.reads(reads1,genes)
consol.map1 <- consolidate.reads(map1)

#reads2 <- load.map.file("test/L1_2394_WT_FRUC_T1.map")
map2 <- map.reads(reads2,genes)
consol.map2 <- consolidate.reads(map2)

find.mutants <- function(maps,n.pick){
  counts <- lapply(maps,function(x) x$reads)
  genes <- lapply(maps,function(x) x$gene)
  total <- sapply(counts,sum)
  total.genes <- length(unique(c(genes,recursive=T)))
  freq <- counts
  for (i in 1:length(freq))
    freq[[i]] <- freq[[i]] / total[i]
  n.libraries <- length(freq)
  counts.sep <- numeric(length(n.pick))
  counts.pooled <- counts.sep
  for (i in 1:length(n.pick)) {
    found.sep <- c()
    found.pooled <- c()
    for (j in 1:n.libraries) {
      found.sep <- unique(c(found.sep, genes[[j]][n.pick[i]*freq[[j]] >= 1]))
      found.pooled <- unique(c(found.pooled, genes[[j]][n.pick[i]*freq[[j]] >= n.libraries]))
    }
    counts.sep[i] <- length(found.sep)
    counts.pooled[i] <- length(found.pooled)
  }
  plot(n.pick/1000,counts.sep/total.genes * 100,
       type="l",
       xlab="# colonies picked (thousands)",
       ylab=paste("% mutants found (", total.genes, " total)", sep=""),
       xaxt="n")
  lines(n.pick/1000,counts.pooled/total.genes * 100,
        col="red")
  axis(1, at = seq(1, 100, by = 10))
}

find.mutants(maps=list(consol.map1,consol.map2), n.pick=seq(1000,100000,1000))

#counts <- consol.map$reads
#total <- sum(counts)
#total.genes <- dim(consol.map)[[1]]
#freq <- counts / total

#n.pick <- seq(1000,100000,1000)
#get.n.found <- function(n) sum(freq*n >= 1)

#n.found <- sapply(n.pick,get.n.found)
#plot(n.pick/1000,n.found/total.genes * 100,
#     type="l",
#     xlab="# colonies picked (thousands)",
#     ylab=paste("% mutants found (", total.genes, " total)", sep=""),
#     xaxt="n")
#axis(1, at = seq(1, 100, by = 10))