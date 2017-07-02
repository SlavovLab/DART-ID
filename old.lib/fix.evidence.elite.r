#+---------------------------------------------------------------------------------------------------------------------+
#| 1.a. FIX.EVIDENCE: (v. 1)                                                                                           |
#|                                                                                                                     |
#| This function fixes the evidence.txt as follows:                                                                    |
#| If we have similar peptides, this function takes the sequence i and:                                                |
#| Sets PEP as min(PEP(i)),                                                                                            |
#| and RT as median(Retention.times(i)).                                                                               |
#| Also evidence.txt will be shrinked and filtered to PEP<0.05 to occupy less memory.                                  |
#| The output ev.ed.txt will be used for making the library and selecting the reference experiment.                    |
#|                                                                                                                     |
#| Example: fix.evidence(path.in = 'C:/../evidence.txt',                                                               |
#|                      path.out = 'C:/../ev.ed.txt')                                                                  |
#|                                                                                                                     |
#+---------------------------------------------------------------------------------------------------------------------+

fix.evidence.elite <- function(path.in, path.out){
  # load data ---------------------------------------------------------------
  
  # load tab-delimited files, with only the columns specified below
  ev <- read_tsv(path.in,
                 col_types=cols_only(
                   Sequence=col_character(),
                   `Raw file`=col_character(),
                   `Retention time`=col_number(),
                   PEP=col_number()
                 ))
  # change column names
  colnames(ev) <- c('Sequence','Raw.file','Retention.time','PEP')
  # only take peptides with PEP < 0.05
  ev <- ev[ev$PEP<0.05,]
  
  # get the list of all experiments (raw files) in the input file
  exps <- unique(ev$Raw.file)
  
  # exclude some experiments here??
  exps <- exps[grep('[0-9]{6}?A', exps)]
  
  # create output table
  all.exps <- data.frame(Raw.file=character(), 
                         Sequence=character(), 
                         PEP=numeric(), 
                         Retention.time=numeric())
  
  # Consolidate Peptides ----------------------------------------------------
  
  # loop thru experiments
  counter <- 0
  for (i in exps) {
    # get current experiment
    exp <- subset(ev, ev$Raw.file==i)
    
    counter <- counter + 1
    cat('\r', 'Processing ', counter, '/', length(exps), ' ', i,  '                                          ')
    flush.console()
    
    # count # of "scans"/IDs for each peptide
    # some of these will come from consolidated PSMs, and will not be further consolidated because
    # of different charges, i.e., +1, +2
    pep.frqs <- as.data.frame(table(exp$Sequence))
    
    # get peptides where we have more than one "scan"
    red.peps.list <- as.character(pep.frqs[pep.frqs$Freq>1, 'Var1'])
    # for these peptides, take the smallest PEP and the median RT
    red.peps <- data.frame(Raw.file=as.character(rep(i, length(red.peps.list))),
                           Sequence=as.character(rep(NA, length(red.peps.list))),
                           PEP=as.numeric(rep(NA, length(red.peps.list))),
                           Retention.time=as.numeric(rep(NA, length(red.peps.list))),
                           stringsAsFactors=FALSE)
    if(length(red.peps.list) > 0) {
      for (j in 1:length(red.peps.list)) {
        red.peps$Sequence[j] <- red.peps.list[j]
        red.peps$PEP[j] <- min(unlist(exp[exp$Sequence==red.peps.list[j], 'PEP']))
        red.peps$Retention.time[j] <- median(unlist(exp[exp$Sequence==red.peps.list[j], 'Retention.time']))
      }
    }
    
    # consolidate list so that all peptides have a unique PEP and RT
    exp <- rbind(red.peps, exp[!(exp$Sequence %in% red.peps.list), ])
    # add to output table
    all.exps <- rbind(all.exps, exp)
  }
  
  
  # Write Output ------------------------------------------------------------
  
  # scrub row names
  row.names(all.exps) <- NULL
  # write to output file - ev.ed.txt
  write.table(all.exps, path.out, sep = '\t', row.names = FALSE, quote = FALSE)
}









