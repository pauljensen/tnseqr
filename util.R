
get_grouping <- function() {
  getOption("tnseqr_groupings")
}

set_grouping <- function(grouping) {
  options(tnseqr_groupings=grouping)
}
