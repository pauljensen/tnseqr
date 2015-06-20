
Rcpp::sourceCpp('hash.cpp')
build_hashes(c("a","b","c"))
add_keys("a", c("x", "y", "z"))
add_keys("b", c("x"))
get_values("a", c("x", "y", "z"))
get_values("b", c("x", "y"))
write_hash_to_file("a", "test.fasta")
erase_hashes(c("a","b","c"))
build_hashes(c("a","b","c"))
get_values("a", c("x", "y", "z"))

