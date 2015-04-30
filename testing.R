
require(stringr)
require(dplyr)
require(tidyr)
require(GenomicRanges)
require(GENE.E)
require(parallel)

# turn off multicore computations by default
options(mc.cores=1)

source("preprocess.R")
source("compile_reads.R")
source("parse_genbank.R")
source("homology.R")
source("plotting.R")
source("sampling.R")


