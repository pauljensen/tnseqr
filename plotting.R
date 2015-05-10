
send_to_genee <- function(tnseq, homologs=NULL) {
  tnseq$fitness$experiment <- paste(tnseq$fitness$condition,
                                    tnseq$fitness$strain)
  exper_annot <- tnseq$fitness %>% select(strain, condition, experiment) %>%
                    group_by(strain, condition, experiment) %>% summarize()
  
  if (!is.null(homologs)) {
    mapping <- make_mapping(homologs)
    tnseq$fitness$locus <- mapping[tnseq$fitness$locus]
  }
  fit <- select(tnseq$fitness, experiment, fitness, locus)
  fit <- fit[!is.na(fit$locus),]
  mat <- spread(fit, experiment, fitness)
  rownames(mat) <- mat$locus
  mat$locus <- NULL
  to.genee(as.matrix(mat))
  return(fit)
}

plot_insertion_density <- function(tnseq, max_insertions=20) {
  data <- tnseq$fitness %>% 
    filter(!is.na(locus), insertions < max_insertions)
  (ggplot(aes(x=insertions), data=data) + 
     geom_density(alpha=0.5) + 
     theme_bw() + 
     facet_grid(strain ~ library)) %>%
    print()
}

plot_condition_strain_matrix <- function(tnseq, variable="fitness") {
  fit <- tnseq$fitness %>% group_by(strain, condition) %>% mutate(xid=1:n())
  (ggplot(aes_string(x="xid", y=variable, color="library"), data=fit) +
     facet_grid(strain ~ condition) + 
     geom_point(alpha=1/20) + 
     theme_bw()) %>%
    print()
}