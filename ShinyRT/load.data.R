load.data <- function(file.handle) {
  library(readr)
  
  ev <- read_tsv(file.handle)
  
  return(ev)
}
