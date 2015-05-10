
load_timseq_file <- function(filename) {
  read.csv(file=filename, stringsAsFactors=F)
}

load_timseq_files <- function(dir, libraries, strain) {
  files <- paste0(dir, "/", libraries, "_", strain, ".csv")
  loaded <- lapply(files, load_timseq_file)
  for (i in seq_along(libraries)) {
    loaded[[i]]$library <- libraries[i]
  }
  rbind_all(loaded)
}

load_timseq_files(dir="~/Dropbox/bc/tnseqr/test/Glucose - T4vsD39vs19F - March 2014/Individual insertions/2394 vs 326-D39-Glucose-IndividInserts",
                  libraries=c("L1", "L2", "L3", "L4", "L5"),
                  strain="326_Glucose") -> data

# load_timseq_files(dir="~/Dropbox/bc/tnseqr/test/Glucose - T4vsD39vs19F - March 2014/Individual insertions/2394 vs 19F -Glucose-IndividInserts",
#                   libraries=c("L1", "L2", "L3", "L4", "L5", "L6"),
#                   strain="T4_Gluc") -> data

timseq_fitness <- data %>% 
  filter(nchar(gene) > 0, !is.na(W)) %>% 
  group_by(gene, library) %>% 
  summarise(fitness=mean(W), 
            fitness_sd=sd(W),
            n_insertions=n())

#t %<>% calculate_fitness()
tnseq_fitness <- t$fitness %>% filter(strain=="D39", condition=="SDMM") %>%
  filter(!is.na(fitness), is.finite(fitness)) %>%
  group_by(locus, library) %>% summarise(fitness_sd=sd(fitness),
                                fitness=mean(fitness),
                                n_insertions=sum(insertions))

joined_d39 <- inner_join(timseq_fitness, tnseq_fitness, by=c("gene"="locus"))


qplot(x=fitness.x, y=fitness.y, data=joined_t4) %>% print()
qplot(x=n_insertions.x, y=n_insertions.y, data=joined_t4) %>% print()
qplot(x=fitness_sd.x, y=fitness_sd.y, data=joined_t4) %>% print()

