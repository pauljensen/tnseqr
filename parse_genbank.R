
load_genbank_file <- function(filename,key="gene") {
  lines <- readLines(filename)
  is_entry <- str_detect(lines, "^     \\w")
  is_keyed <- str_detect(lines, paste0("^     ", key))
  starts <- which(is_keyed)
  n_groups <- length(starts)
  stops <- integer(n_groups)
  all_stops <- which(c(is_entry, T))
  for (i in 1:n_groups) {
    stops[i] <- (all_stops[all_stops > starts[i]])[1] - 1
  }
  groups <- lapply(1:n_groups, function(i) lines[starts[i]:stops[i]])
  
  parse_group <- function(group) {
    range <- as.integer(str_match(group[1], "(\\d+)..(\\d+)")[1,c(2,3)])
    get_quoted <- function(name) {
      matches <- str_match(group, paste0("/", name, "=\"(.+)\""))
      hits <- matches[!is.na(matches[,2]),2]
      if (length(hits) == 0) {
        return(NA)
      } else if (length(hits) > 1) {
        return(hits[1])
      } else {
        return(hits)
      }
    }
    return(data.frame(start=range[1], stop=range[2],
                      locus=get_quoted("locus_tag"),
                      gene=get_quoted("gene"),
                      stringsAsFactors=F))
  }
  
  df <- NULL
  do.call(rbind, lapply(groups, parse_group))
}

