
require(GENE.E)
require(tidyr)

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

#to.genee(plot_fieldvalue(cfit, "condition", "CDM", "strain"))
