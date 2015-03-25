
plot_fieldvalue <- function(fitness, pull, value, cols) {
  fitness %<>% group_by(strain, condition, locus) %>% 
      summarise(fitness=mean(fitness), homolog=homolog[1])
  d <- fitness[fitness[[pull]] == value, c("homolog", "fitness", cols)]
  wider <- as.data.frame(spread_(d, cols, "fitness"))
  rownames(wider) <- wider$homolog
  wider$homolog <- NULL
  
  #heatmap.2(as.matrix(wider) - 1, dendrogram="none", trace="none", Colv=F, Rowv=F,
  #          col=colorRampPalette(c(rgb(1,0,0), rgb(1,1,1), rgb(0,0,1))),
  #          scale="none", key=F
  #)
  return(as.matrix(wider))
}

w <- plot_fieldvalue(fit, "condition", "CDM", "strain")

ggplot(cfit %>% filter(condition=="CDM"), 
       aes(x=homolog, y=strain, fill=fitness)) + 
  geom_tile() + theme_bw() -> p
  #scale_colour_gradientn(colours=brewer.pal(10,"RdBu")) -> p
#print(p)


to.genee(plot_fieldvalue(cfit, "condition", "CDM", "strain"))
