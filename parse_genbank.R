
library(plyr)

get.word.starts <- function(line) {
  gregexpr('[^ ]+',line,perl=TRUE)[[1]]
}

parse.field.line <- function(line) {
  match <- gregexpr('[^ ]+',line,perl=TRUE)[[1]]
  list(
    name = substr(line, match[1], match[1]+attr(match,'match.length')[1]-1),
    name.start = match[1],
    value = substring(line, match[2]),
    value.start = match[2],
    indent = match[2] - 1
  )
}

parse.block.indented <- function(lines) {
  index <- 0
  n.lines <- length(lines)
  
  next.line <- function() {
    index <<- index + 1
    lines[index]
  }
  #current.line <- function() lines[index]
  backstep.line <- function() index <<- index - 1
  peek.line <- function() lines[index+1]
  is.another.line <- function() index < n.lines
  
  finish.field <- function(value.start,value) {
    while (is.another.line() && get.word.starts(peek.line())[1] == value.start) {
      line <- next.line()
      starts <- get.word.starts(line)
      value <- c(value, substring(line,value.start))
    }
    value
  }
  
  get.subfield <- function() {
    field.match <- parse.field.line(next.line())
    value <- field.match$value
    value <- finish.field(field.match$value.start, field.match$value)
    list(
      name=field.match$name,
      value=value,
      fields=list()
    )
  }
  
  get.next.field <- function() {
    # assumes there is a next line
    field.match <- parse.field.line(next.line())
    field <- list(name=field.match$name,
                  value=field.match$value,
                  fields=list())
    while (is.another.line() && get.word.starts(peek.line())[1] > field.match$name.start) {
      line <- next.line()
      starts <- get.word.starts(line)
      if (starts[1] == field.match$value.start) {
        # continuation of previous value
        field$value <- c(field$value, substring(line,starts[1]))
      } else {
        # subfield
        backstep.line()
        field$fields <- c(field$fields, list(get.subfield()))
      }
    }
    field
  }
  
  fields <- list()
  while (is.another.line()) {
    fields <- c(fields, list(get.next.field()))
  }
  fields
}

fields.to.data.frame <- function(fields) {
  parse.qualifier <- function(qual) {
    matches <- regmatches(qual, regexec("/(.+)=(.+)",qual))[[1]][-1]
    if (substr(matches[2],1,1) == '"')
      matches[2] <- substr(matches[2],2,nchar(matches[2])-1)
    matches
  }
  build.row <- function(field) {
    field$value[1] <- paste("/range=", field$value[1], sep="")
    qual.list <- lapply(as.list(grep('^/.*=',field$value,perl=T,value=T)), parse.qualifier)
    quals <- vapply(qual.list, function(x) x[2], "")
    names(quals) <- vapply(qual.list, function(x) x[1], "")
    as.data.frame(as.list(quals),stringsAsFactors=FALSE)
  }
  lapply(fields, build.row)
  do.call(rbind.fill, lapply(fields, build.row))
}

add.start.stop <- function(genes) {
  start.stop <- lapply(as.list(genes$range), function(x) regmatches(x, regexec("(\\d+)..(\\d+)",x))[[1]][-1])
  genes$start <- as.integer(vapply(start.stop, function(x) x[1], ""))
  genes$stop <- as.integer(vapply(start.stop, function(x) x[2], ""))
  genes
}

parse.genbank.lines <- function(lines) {
  # returns a data frame of the genes
  fields <- parse.block.indented(lines)
  names(fields) <- vapply(fields, function(x) x$name, "")
  genes <- fields.to.data.frame(Filter(function(x) x$name == "gene", fields$FEATURES$fields))
  add.start.stop(genes)
}

parse.genbank.file <- function(filename) parse.genbank.lines(readLines(filename))




genes <- parse.genbank.file("test/tigr4.gb")



