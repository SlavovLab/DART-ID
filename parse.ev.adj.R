parse.ev.adj <- function(file.path, type='STAN') {
  library(readr)
  
  if(is.character(file.path)) {
    ev <- read_tsv(file.path)
  } else {
    ev <- file.path
  }
  
  if(type == 'STAN') {
    ev <- ev[, c('Sequence', 'Proteins', 'Raw file', 
                 'Retention time','PEP', 'rt.minus', 'rt.plus',
                 'muijs', 'sigmas', 'PEP.new', 
                 "Reporter intensity corrected 0", "Reporter intensity corrected 1", 
                 "Reporter intensity corrected 2", "Reporter intensity corrected 3",
                 "Reporter intensity corrected 4", "Reporter intensity corrected 5", 
                 "Reporter intensity corrected 6", "Reporter intensity corrected 7", 
                 "Reporter intensity corrected 8", "Reporter intensity corrected 9")]
  } else if(type == 'Ali') {
    ev <- ev[, c('Sequence', 'Leading.razor.protein', 'Raw.file',
                 'Retention.time', 'PEP', 'RT.lib', 'RT.corrected', 
                 'dRT', 'PEP.new',
                 "Reporter.intensity.corrected.0", "Reporter.intensity.corrected.1", 
                 "Reporter.intensity.corrected.2", "Reporter.intensity.corrected.3",
                 "Reporter.intensity.corrected.4", "Reporter.intensity.corrected.5", 
                 "Reporter.intensity.corrected.6", "Reporter.intensity.corrected.7", 
                 "Reporter.intensity.corrected.8", "Reporter.intensity.corrected.9")]
  }
  
  # PEP > 1 defaults to PEP = 1
  ev$PEP[ev$PEP > 1] <- 1
  
  # PEP.updated - PEP.new except when its NA, and then defaults to PEP
  ev$PEP.updated <- ev$PEP.new
  ev$PEP.updated[is.na(ev$PEP.new)] <- ev$PEP[is.na(ev$PEP.new)]
  # PEP == 0 defaults to PEP = min(PEP)
  ev$PEP.updated[ev$PEP.updated <= 0] <- min(ev$PEP.updated)
  
  # updated RT
  ev$RT.new <- ev$muijs
  ev$RT.new[is.na(ev$RT.new)] <- ev$`Retention time`[is.na(ev$RT.new)]
  
  # delta PEP and RT
  ev$dPEP <- ev$PEP.new - ev$PEP
  if(type=='STAN') {
    ev$dRT <- ev$muijs - ev$`Retention time`
  }
  
  # calculate q values
  # as cumulative average (running average) of the PEP
  # (see PEP and FDR - two sides of the same coin)
  ev$q <- (cumsum(sort(ev$PEP.updated)) / seq(1,nrow(ev)))[order(ev$PEP.updated)]
  ev$q.old <- (cumsum(sort(ev$PEP)) / seq(1,nrow(ev)))[order(ev$PEP)]
  
  ev$ID <- seq(1,nrow(ev))
  
  return(ev)
}