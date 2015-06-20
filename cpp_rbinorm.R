
require(Rcpp)

cppFunction('NumericVector sampleReads(NumericVector x, double frac) {
  int n = x.size();
  NumericVector out(n);

  for(int i = 0; i < n; ++i) {
    out[i] = rbinom(1,x[i],frac)[0];
  }
  return out;
}')

sampleReadsR <- function(reads, frac) {
  vapply(reads, function(x) rbinom(1,x,frac), 0)
}

sampleReadsFor <- function(reads, frac) {
  for (i in seq_along(reads))
    reads[i] <- rbinom(1,reads[i],frac)
  reads
}

reads <- as.integer(runif(10000)*10000)
frac <- 0.1
microbenchmark(
  sampleReads(reads, frac),
  sampleReadsR(reads, frac),
  sampleReadsFor(reads, frac)) -> results
print(results)
