# clean up evidence files so that they can be converted to mzTab by MassIVE

library(tidyverse)

folder_path <- '/gd/bayesian_RT/MQ_evidence'
files <- list.files(folder_path)

for (f in files) {
  ev <- read_tsv(file.path(folder_path, f))
  
  # remove any match-between-runs hits
  ev <- ev %>%
    filter(!is.na(PEP))
  
  ev <- ev %>%
    # remove leading and trailing '_' from modified sequence
    mutate(`Modified sequence`=gsub('_', '', `Modified sequence`)) %>%
    # ceil PEPs to 1
    mutate_at('PEP', funs(ifelse(. > 1, 1, .))) %>%
    # calculate q-values
    mutate(qval=(cumsum(PEP[order(PEP)]) / 
                   seq(1, nrow(ev)))[order(order(PEP))])
  
  # write to file
  write_tsv(ev, file.path(folder_path, paste0('updated_', f)))
}