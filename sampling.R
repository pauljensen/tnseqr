
calc_insertion_stats <- function(tnseq) {
  inserts_before <- tnseq$insertions %>% group_by(strain, condition) %>% summarize(insertions_before=n())
  filtered <- tnseq %>% filter_insertions()
  inserts_after <- filtered$insertions %>% group_by(strain, condition) %>% summarize(insertions_after=n())
  counts <- inner_join(inserts_before, inserts_after, by=c("strain", "condition"))
    
  fitness <- filtered %>% calculate_fitness()
  fitsd <- fitness$insertions %>% group_by(strain, condition) %>% summarize(sd_fitness=mean(fitness_sd, na.rm=T),
                                                                            mean_insertions=mean(insertions, na.rm=T),
                                                                            median_insertions=floor(median(insertions, na.rm=T)))
  return(inner_join(counts, fitsd, by=c("strain", "condition")))
}

downsample_insertions <- function(tnseq, fracs=c(0.1, 0.3, 0.5, 0.7, 0.8, 0.9, 0.95)) {
  downsample_aux <- function(frac) {
    cat(frac, "\n")
    t <- tnseq
    binom_reads <- function(cnt) {
      vapply(cnt, function(x) rbinom(1,x,frac), 0)
    }
    t$insertions$reads1 <- binom_reads(t$insertions$reads1)
    t$insertions$reads2 <- binom_reads(t$insertions$reads2)
    
    t$insertions <- t$insertions %>% filter(reads1 > 0 & reads2 > 0)
    
    stats <- calc_insertion_stats(t)
    stats$frac <- frac
    return(stats)
  }
  
  all_stats <- rbind_all(mclapply(fracs, downsample_aux))
  full_stats <- calc_insertion_stats(tnseq)
  full_stats$frac <- 1.0
  all_stats <- bind_rows(all_stats, full_stats)
  
  return(all_stats)
}

#ptm <- proc.time()
#stats <- downsample_insertions(tin, fracs=c(0.1, 0.3, 0.5))
#elapsed <- proc.time() - ptm
#print(elapsed)
