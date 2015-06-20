
calc_insertion_stats <- function(tnseq, grouping=get_grouping()) {
  inserts_before <- tnseq$insertions %>% group_by_(.dots=grouping) %>% summarize(insertions_before=n())
  filtered <- tnseq %>% filter_insertions()
  inserts_after <- filtered$insertions %>% group_by_(.dots=grouping) %>% summarize(insertions_after=n())
  counts <- inner_join(inserts_before, inserts_after, by=c("strain", "condition", "library"))
    
  fitness <- filtered %>% calculate_fitness()
  fitsd <- fitness$fitness %>% group_by_(.dots=grouping) %>% summarize(sd_fitness=mean(fitness_sd, na.rm=T),
                                                             mean_insertions=mean(insertions, na.rm=T),
                                                             median_insertions=floor(median(insertions, na.rm=T)))
  return(inner_join(counts, fitsd, by=grouping))
}

cppFunction('NumericVector sampleReads(NumericVector x, double frac) {
  int n = x.size();
  NumericVector out(n);
  
  for(int i = 0; i < n; ++i) {
    out[i] = rbinom(1,x[i],frac)[0];
  }

  return out;
}')

downsample_insertions <- function(tnseq, fracs=c(0.1, 0.3, 0.5, 0.7, 0.8, 0.9, 0.95)) {
  downsample_aux <- function(frac) {
    cat(frac, "\n")
    t <- tnseq
    
    t$insertions$reads1 <- sampleReads(t$insertions$reads1, frac)
    t$insertions$reads2 <- sampleReads(t$insertions$reads2, frac)
    
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



# ================= insertion dropout ======================

cppFunction('NumericVector sampleFrequencies(NumericVector freq, double total) {
  int n = freq.size();
  NumericVector out(n);
  
  for(int i = 0; i < n; ++i) {
    out[i] = rbinom(1,total,freq[i])[0];
  }

  return out;
}')

plot_insertion_dropout <- function(inserts,reads="reads2",cutoff=17,
                                   fracs=seq(0.1,2.0,0.1),
                                   title_prefix="Library") {
  counts <- inserts[[reads]]
  counts <- counts[counts > 0]
  total <- sum(counts)
  freq <- counts / total
  passing <- numeric(length(fracs))
  print(c(total))
  for (i in seq_along(fracs)) {
    passing[i] <- sum(sampleFrequencies(freq,floor(fracs[i]*total)) > cutoff) / length(freq)
    #cat(i)
  }
  plot(fracs*total,passing,ylim=c(0,1),xlab="reads per sample",
       ylab="fraction sites passing",type="l")
  title(sprintf("%s: %i unique insertions",title_prefix, length(freq)))
  abline(v=total,col="red",lty="dashed")
}

par(mfrow=c(3,2))

for (st in c("316","D39","19F")) {
  for (lib in c("L1","L2","L3","L4","L5","L6")) {
    t$insertions %>% 
      filter(strain==st,library==lib,condition=="CDM") %>% 
      plot_insertion_dropout(title_prefix=paste(st, lib, "CDM"))
  }
}



#ptm <- proc.time()
#stats <- downsample_insertions(tin, fracs=c(0.1, 0.3, 0.5))
#elapsed <- proc.time() - ptm
#print(elapsed)
