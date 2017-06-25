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

fix.evidence <- function(path.in, path.out){
  # read input file as tab-delimited, text file
  ev <- read.delim(path.in, header = TRUE, stringsAsFactors = FALSE)
  # only take peptides with PEP < 0.05
  # only take columns we need
  ev <- ev[ev$PEP<0.05, c('Raw.file', 'Sequence', 'PEP', 'Retention.time')]
  
  # get the list of all experiments (raw files) in the input file
  exps <- unique(ev$Raw.file)
  
  # create output table
  all.exps <- data.frame(Raw.file=character(), Sequence=character(), PEP=numeric(), Retention.time=numeric())
  
  # loop thru experiments
  for (i in exps) {
    # get current experiment
    exp <- subset(ev, ev$Raw.file==i)
    
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
    
    for (j in 1:length(red.peps.list)) {
      red.peps$Sequence[j] <- red.peps.list[j]
      red.peps$PEP[j] <- min(exp[exp$Sequence==red.peps.list[j], 'PEP'])
      red.peps$Retention.time[j] <- median(exp[exp$Sequence==red.peps.list[j], 'Retention.time'])
    }
    
    # consolidate list so that all peptides have a unique PEP and RT
    exp <- rbind(red.peps, exp[!(exp$Sequence %in% red.peps.list), ])
    # add to output table
    all.exps <- rbind(all.exps, exp)
  }
  # scrub row names
  row.names(all.exps) <- NULL
  # write to output file - ev.ed.txt
  write.table(all.exps, path.out, sep = '\t', row.names = FALSE, quote = FALSE)
}









