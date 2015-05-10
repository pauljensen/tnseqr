
require(stringr)
require(dplyr)
require(tidyr)
require(GenomicRanges)
require(GENE.E)
require(parallel)
require(Rcpp)
require(ggplot2)

source("util.R")
source("preprocess.R")
source("compile_reads.R")
source("parse_genbank.R")
#source("homology.R")
source("plotting.R")
source("sampling.R")

# turn off multicore computations by default
options(mc.cores=1)

# pool libraries
set_grouping(c("strain", "library", "condition"))
