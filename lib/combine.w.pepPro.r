
#+---------------------------------------------------------------------------------------------------------------------+
#| 1.b. COMBINE.PEPPRO: (v. 1)                                                                                         |
#|                                                                                                                     |
#| This function combines PeptideProphet processed results with the evidence file and creates ev.ed.pepPro.txt which is|
#| the input to find.ref.exp and make.lib functions. It needs evidence.txt and for peptideProphet processed results it |
#| is assumed that all of the results are in one file, and in TPP output format. This function cleans up the sequences |
#| and removes the numbers and flanking residuals from the sequences of peptides. The file must have the following     |
#| columns with exact same headers:                                                                                    |
#| Raw.file = The exact name of the Raw file of the experiment as it appears in the evidence file. This columns should |
#| be added manually by the user. There is also an option in xinteract to add labels for the experiments but the header|
#| of the column should be exactly as it is stated here.                                                               |
#| peptide = The sequence of peptides as they appear in the output of PepPro processed file. This function cleans the  |
#| sequences using regular expression.                                                                                 |
#| probability = Peptide Prophet probability. Min = 0.8                                                                |
#| retention_time_sec = Retention time in seconds as it is reported by PeptideProphet.                                 |
#| Other columns might also be in the file but will not be used in this function.                                      |
#| The evidence.txt file must also be fed to the function without any alterations.                                     |
#|                                                                                                                     |
#| Example: combine.pepPro(path.in.pepPro = 'C:/../peptideProphet.txt',                                                |
#|                        path.in.evdence = 'C:/../evidence.txt',                                                      |
#|                               path.out = 'C:/../ev.ed.pepPro.txt')                                                  |
#|                                                                                                                     |
#+---------------------------------------------------------------------------------------------------------------------+

combine.w.pepPro <- function(path.in.pepPro, path.in.evidence, path.out){
  all.exp <- read.delim(path.in.pepPro, sep='\t', stringsAsFactors=F, header=T)
  all.exp <- all.exp[, c("Raw.file", "peptide", "probability", "retention_time_sec")]
  all.exp <- all.exp[all.exp$probability>=0.9, ]
  all.exp$peptide <- substr(all.exp$peptide, 3, nchar(all.exp$peptide)-2)
  all.exp$peptide <- gsub("[0-9a-z.]+", "", all.exp$peptide)
  all.exp$peptide <- gsub("\\[\\]", "", all.exp$peptide)
  all.exp$retention_time_sec <- all.exp$retention_time_sec/60
  colnames(all.exp)[2] <- "Sequence"
  colnames(all.exp)[3] <- "PEP"
  colnames(all.exp)[4] <- "Retention.time"
  all.exp[all.exp$PEP>0.999, "PEP"] <- 1e-4
  all.exp[!all.exp$PEP==1e-4, "PEP"] <- 3e-2
  
  ev <- read.delim(path.in.evidence, header = TRUE, stringsAsFactors = FALSE)
  ev <- ev[ev$PEP<0.05, c('Raw.file', 'Sequence', 'PEP', 'Retention.time')]
  
  ev <- rbind(ev, all.exp)
  exps <- unique(ev$Raw.file)
  
  all.exps <- data.frame(Raw.file=character(), Sequence=character(), PEP=numeric(), Retention.time=numeric())
  for (i in exps) {
    exp <- subset(ev, ev$Raw.file==i)
    pep.frqs <- as.data.frame(table(exp$Sequence))
    red.peps.list <- as.character(pep.frqs[pep.frqs$Freq>1, 'Var1'])
    if (length(red.peps.list)==0){
      all.exps <- rbind(all.exps, exp)
      next
    }
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
    
    exp <- rbind(red.peps, exp[!(exp$Sequence %in% red.peps.list), ])
    all.exps <- rbind(all.exps, exp)
  }
  row.names(all.exps) <- NULL
  
  write.table(all.exps, path.out, sep = '\t', row.names = FALSE, quote = FALSE)
}

