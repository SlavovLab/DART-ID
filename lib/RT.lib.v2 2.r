#+---------------------------------------------------------------------------------------------------------------------+
#| Making the library has 4 steps which are accomplished by these five functions:                                      |
#| 1. a.fix.evidence OR b.combine.pepPro : Prepares input of the library maker function.                               |
#| 2. find.ref.exp : Finds the reference experiment.                                                                   |
#| 3. make.lib : Makes the peptide retention time library.                                                             |
#| 4. dRT : Updates the PEPs according to dRTs.                                                                        |
#|                                                                                                                     |
#|                                                                                                                     |
#| All of the inputs of the functions must be tab delimited files. The outputs are all tab delimited as well.          |
#|                                                                                                                     |
#| -Ali Banijamali                                                                                                     |
#+---------------------------------------------------------------------------------------------------------------------+
#
#
#
#
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







#+---------------------------------------------------------------------------------------------------------------------+
#| 2. FIND.REF.EXP: (v. 1)                                                                                             |
#|                                                                                                                     |
#| This function finds the best experiment as ref.exp. The ref.exp will be printed after the function finishes but a   |
#| txt file will also be saved which has information about how all of the experiments fit. We have to manually give    |
#| the name of the reference experiment to the make.lib function.                                                      |
#| This function takes the fixed evidence file prepared by fix.evidence function as input.                             |
#| The output contains These info for each experiment:                                                                 |
#| c.fit: Number of exps with less than 10 shared peptides between ref and new exp.                                    |
#| c.med.res: Number of exps with >2 min median residuals of the fit.                                                  |
#| c.mean.res: Number of exps with >2 min mean residuals of the fit.                                                   |
#| score: The experiment with lowest score (c.fit+c.med.res+c.mean.res) is ref.exp.                                    |
#|                                                                                                                     |
#| This function needs 'MASS' package for the robust regression and it is assumed that it is installed and required.   |
#|                                                                                                                     |
#| Example: find.ref.exp(path.in = 'C:/../ev.ed.txt',                                                                  |
#|                      path.out = 'C:/../ref.picks.txt')                                                              |
#|                                                                                                                     |
#+---------------------------------------------------------------------------------------------------------------------+

find.ref.exp <- function(path.in, path.out){
  #install.packages('MASS')
  #library('MASS')
  
  # import fixed evidence file  
  ev.ed <- read.delim(path.in, header = TRUE, stringsAsFactors = FALSE)
  
  # only include experiments that have peptides with PEP < 5e-4
  inc.list <- unique(subset(ev.ed, ev.ed$PEP<5e-4)$Raw.file)
  ev.ed <- ev.ed[ev.ed$Raw.file %in% inc.list, ]
  
  #SELECTING THE REFERENCE EXPERIMENT:
  # create frame for candidates
  ref.pick.candidates <- data.frame(experiment=unique(ev.ed$Raw.file),
                                    c.fit=rep(NA, length(unique(ev.ed$Raw.file))),
                                    c.med.res=rep(NA, length(unique(ev.ed$Raw.file))),
                                    c.mean.res=rep(NA, length(unique(ev.ed$Raw.file))))
  
  # loop thru candidates
  exps.all <- unique(ev.ed$Raw.file)
  for (ii in exps.all) {
    C10 <- 0 
    
    # load reference experiment
    reference.exp.l <- subset(ev.ed, ev.ed$PEP<5e-4 & ev.ed$Raw.file==ii)
    # sort it alphabetically by peptide sequence
    reference.exp.l <- reference.exp.l[with(reference.exp.l, order(Sequence)), ]
    
    # scrub row names
    row.names(reference.exp.l) <- NULL
    # get experiments that aren't the one we are analyzing
    exps <- exps.all[!exps.all %in% ii]
    
    # create coefficients frame
    shift.coeffs <- data.frame(experiment=unique(ev.ed$Raw.file),
                               intercept=rep(NA, length(unique(ev.ed$Raw.file))),
                               coeff=rep(NA, length(unique(ev.ed$Raw.file))))
    shift.coeffs[shift.coeffs$experiment==ii, 2:3] <- 1
    
    # compare this experiment to all others
    for (i in exps) {
      # load comparison
      new.exp <- subset(ev.ed, ev.ed$Raw.file==i)
      row.names(new.exp) <- NULL
      
      # only use comparison with peptides with PEP < 5e-4, and with peptides in common with
      # the experiment we are currently analyzing
      new.exp.f <- subset(ev.ed, ev.ed$Raw.file==i & ev.ed$PEP<5e-4)
      new.exp.f <- new.exp.f[new.exp.f$Sequence %in% reference.exp.l$Sequence, ]
      new.exp.f <- new.exp.f[with(new.exp.f, order(Sequence)), ]
      row.names(new.exp.f) <- NULL
      new.exp.f <- droplevels(new.exp.f)
      
      # get data of shared peptides from current experiment
      ref.shared <- reference.exp.l[reference.exp.l$Sequence %in% new.exp.f$Sequence, ]
      row.names(ref.shared) <- NULL
      ref.shared <- droplevels(ref.shared)
      
      # if more than 10 shared peptides, then increment C10 by 1
      if (nrow(ref.shared)<10){ 
        C10 <- C10+1
        next
      }
      
      # fit robust linear model  
      fit <- rlm(ref.shared$Retention.time ~ new.exp.f$Retention.time, psi=psi.bisquare, maxit=100)
      shift.coeffs[shift.coeffs$experiment==i, 2:3] <- fit$coefficients
      
      ## calculate residuals from the RLM
      # adjust retention times based on fit coefficients (slope and intercept)
      fit.residuals <- as.matrix(data.frame(rep(1, nrow(new.exp.f)), new.exp.f$Retention.time)) %*% t(as.matrix(shift.coeffs[shift.coeffs$experiment==i, 2:3]))
      # re-match residuals to peptide sequences
      fit.residuals <- data.frame(peptides=new.exp.f$Sequence, shifted.RTs=as.numeric(fit.residuals))
      fit.residuals$ref.RTs <- ref.shared$Retention.time[match(fit.residuals$peptides, ref.shared$Sequence)]
      # calculate residuals (model - measured)
      fit.residuals$residuals <- abs(fit.residuals$shifted.RTs-fit.residuals$ref.RTs)
      
      # calculate mean and median residuals
      shift.coeffs[shift.coeffs$experiment==i, "mean.rt.res"] <- as.numeric(mean(fit.residuals$residuals))
      shift.coeffs[shift.coeffs$experiment==i, "median.rt.res"] <- as.numeric(median(fit.residuals$residuals))
    }
    
    # consolidate metrics for this experiment
    ref.pick.candidates[ref.pick.candidates$experiment==ii, "c.fit"] <- C10
    ref.pick.candidates[ref.pick.candidates$experiment==ii, "c.med.res"] <- nrow(subset(shift.coeffs, shift.coeffs$median.rt.res>2))
    ref.pick.candidates[ref.pick.candidates$experiment==ii, "c.mean.res"] <- nrow(subset(shift.coeffs, shift.coeffs$mean.rt.res>2))
  }
  
  # The experiment with lowest score is ref.exp
  ref.pick.candidates$score <- ref.pick.candidates$c.fit+ref.pick.candidates$c.med.res+ref.pick.candidates$c.mean.res
  ref.pick <- as.character(ref.pick.candidates[ref.pick.candidates$score==min(ref.pick.candidates$score), "experiment"])
  
  print(ref.pick)
  
  # write to output
  write.table(ref.pick.candidates, path.out, sep = '\t', row.names = FALSE, quote = FALSE)
}







#+---------------------------------------------------------------------------------------------------------------------+
#| 3. MAKE.LIB: (v. 2)                                                                                                 |
#|                                                                                                                     |
#| This function makes the retention time library. It takes the following inputs in it's command:                      |
#| 1. The path to the edited evidence made by fix.evidence function.                                                   |
#| 2. The full name of the reference experiment as reported by find.ref.exp.                                           |
#| 3. The path for saving the library & Retention.time shifts info for experiments. We will need both of these files   |
#| for calculating the retention times of the peptide for each experiment.                                             |
#| The output is the Retention time library and the shifting coefficients for the experiments in the evidence file.    |
#| In the library, for each peptide mean, median, standard deviation and count of RTs is calculated.                   |
#|                                                                                                                     |
#| This function needs 'MASS' package for the robust regression and it is assumed that it is installed and required.   |
#|                                                                                                                     |
#| Example: make.lib(path.in = 'C:/../ev.ed.txt',                                                                      |
#|              path.out.lib = 'C:/../RT.lib.txt',                                                                     |
#|           path.out.coeffs = 'C:/../shift.coeffs.txt',                                                               |
#|                   ref.exp = '160623AEL160624PMRSAM00092_NC126')                                                     |
#| v. 2: Changes in how the function removes outliers.                                                                 |
#+---------------------------------------------------------------------------------------------------------------------+

make.lib <- function(path.in, path.out.lib, path.out.coeffs, ref.exp){
  #install.packages('MASS')
  #library('MASS')
  
  # load inputs
  ev.ed <- read.delim(path.in, header = TRUE, stringsAsFactors = FALSE)
  ref.pick <- ref.exp
  
  # load list of experiments
  exps <- unique(ev.ed$Raw.file)
  
  # load reference experiment data, and sort by peptide sequence
  reference.exp.l <- subset(ev.ed, ev.ed$PEP<5e-4 & ev.ed$Raw.file==ref.pick)
  reference.exp.l <- reference.exp.l[with(reference.exp.l, order(Sequence)), ]
  row.names(reference.exp.l) <- NULL
  
  # remove reference experiment from master list
  exps <- exps[!exps %in% ref.pick]
  
  # BUILDING THE LIBRARY:
  # initialize RT library data frame
  rt.lib <- data.frame(peptides=as.character(sort(unique(ev.ed$Sequence))),
                       matrix(NA, nrow=length(unique(ev.ed$Sequence)), ncol=length(unique(ev.ed$Raw.file))),
                       stringsAsFactors=FALSE)
  
  # load reference experiment but with less stringent threshold for PEP (0.05)
  reference.exp.h <- subset(ev.ed, ev.ed$PEP<0.05 & ev.ed$Raw.file==ref.pick)
  
  # load retention times in from reference experiment
  rt.lib$X1 <- reference.exp.h$Retention.time[match(rt.lib$peptides, reference.exp.h$Sequence)]
  
  # initialize shift data frame
  shift.coeffs <- data.frame(experiment=unique(ev.ed$Raw.file),
                             intercept=rep(NA, length(unique(ev.ed$Raw.file))),
                             coeff=rep(NA, length(unique(ev.ed$Raw.file))))
  # normalize to the reference experiment
  shift.coeffs[shift.coeffs$experiment==ref.pick, 2:3] <- 1
  
  # counters
  C10 <- 0
  C.res <- 0
  counter <- 2
  
  # for every experiment that is not the reference:
  for (i in exps) {
    counter <- counter+1
    
    # grab the experiment
    new.exp <- subset(ev.ed, ev.ed$Raw.file==i)
    row.names(new.exp) <- NULL
    
    # only take confidently identified peptides, sharing sequences w/ reference experiment, then sort
    new.exp.f <- subset(ev.ed, ev.ed$Raw.file==i & ev.ed$PEP<5e-4)
    new.exp.f <- new.exp.f[new.exp.f$Sequence %in% reference.exp.l$Sequence, ]
    new.exp.f <- new.exp.f[with(new.exp.f, order(Sequence)), ]
    row.names(new.exp.f) <- NULL
    new.exp.f <- droplevels(new.exp.f)
    
    # get data for these shared & confidently identified peptides
    ref.shared <- reference.exp.l[reference.exp.l$Sequence %in% new.exp.f$Sequence, ]
    row.names(ref.shared) <- NULL
    ref.shared <- droplevels(ref.shared)
    
    # iterate counter if less than 10 shared peptides
    # then continue loop to next experiment
    if (nrow(ref.shared)<10){ 
      C10 <- C10+1
      rt.lib[, counter] <- NA
      next
    }
    
    # build linear model and store coefficients
    fit <- rlm(ref.shared$Retention.time ~ new.exp.f$Retention.time, psi=psi.bisquare, maxit=100)
    shift.coeffs[shift.coeffs$experiment==i, 2:3] <- fit$coefficients
    
    fit.residuals <- as.matrix(data.frame(rep(1, nrow(new.exp.f)), new.exp.f$Retention.time)) %*% t(as.matrix(shift.coeffs[shift.coeffs$experiment==i, 2:3]))
    fit.residuals <- data.frame(peptides=new.exp.f$Sequence, shifted.RTs=as.numeric(fit.residuals))
    fit.residuals$ref.RTs <- ref.shared$Retention.time[match(fit.residuals$peptides, ref.shared$Sequence)]
    fit.residuals$residuals <- abs(fit.residuals$shifted.RTs-fit.residuals$ref.RTs)
    shift.coeffs[shift.coeffs$experiment==i, "mean.rt.res"] <- as.numeric(mean(fit.residuals$residuals))
    shift.coeffs[shift.coeffs$experiment==i, "median.rt.res"] <- as.numeric(median(fit.residuals$residuals))
    shift.coeffs[shift.coeffs$experiment==ref.pick, "mean.rt.res"] <- 0
    shift.coeffs[shift.coeffs$experiment==ref.pick, "median.rt.res"] <- 0
    
    edited.rts <- as.matrix(data.frame(rep(1, nrow(new.exp)), new.exp$Retention.time)) %*% t(as.matrix(shift.coeffs[shift.coeffs$experiment==i, 2:3]))
    edited.rts <- data.frame(peptides=new.exp$Sequence, shifted.RTs=as.numeric(edited.rts))
    #edited.rts$initial.RTs <- new.exp$Retention.time
    #edited.rts$residuals <- abs(edited.rts$shifted.RTs-edited.rts$initial.RTs)
    
    if (shift.coeffs[shift.coeffs$experiment==i, "median.rt.res"]>2){ 
      C.res <- C.res+1
      rt.lib[, counter] <- NA
      next
    }
    
    rt.lib[, counter] <- edited.rts$shifted.RTs[match(rt.lib$peptides, edited.rts$peptides)]
  }
  
  # LIBRARY POST PROCESSING:
  
  rt.lib[rt.lib<0] <- NA
  
  rt.lib$RT.count <- apply(rt.lib[,2:(ncol(rt.lib))], 1, function(x) sum(!is.na(x)))
  rt.lib <- rt.lib[!rt.lib$RT.count==0,]
  row.names(rt.lib) <- NULL
  
  rt.lib$rt.median <- apply(rt.lib[,2:(length(unique(ev.ed$Raw.file))+1)], 1, median, na.rm=T)
  
  rt.lib$rt.mean <- apply(rt.lib[,2:(length(unique(ev.ed$Raw.file))+1)], 1, mean, na.rm=T)
  
  rt.lib$sd.pep <- apply(rt.lib[,2:(length(unique(ev.ed$Raw.file))+1)], 1, sd, na.rm=T)
  rt.lib[is.na(rt.lib$sd.p), 'sd.pep'] <- 0
  
  write.table(rt.lib, path.out.lib, sep = '\t', row.names = FALSE, quote = FALSE)
  write.table(shift.coeffs, path.out.coeffs, sep = '\t', row.names = FALSE, quote = FALSE)
}







#+---------------------------------------------------------------------------------------------------------------------+
#| 4. dRT (v. 2)                                                                                                       |
#|                                                                                                                     |
#| This function shifts the retention times of experiments, calculates dRTs and updates PEPs                           |
#| As input, it takes the path to library file, shift coefficients file, evidence.txt file from Max Quant and the      |
#| output will be the evidence file with updated PEPs and containing experiments and peptides that were present in the |
#| library. Contaminants are removed during the process.                                                               |
#|                                                                                                                     |
#| Example: dRT(path.in.lib = 'C:/../RT.lib.txt',                                                                      |
#|           path.in.coeffs = 'C:/../shift.coeffs.txt',                                                                |
#|         path.in.evidence = 'C:/../evidence.txt',                                                                    |
#|     path.out.PEP.updated = 'C:/../evidence+dRT.txt')                                                                |
#|                                                                                                                     |
#|                                                                                                                     |
#|For large dRTs PEP.new sometimes gets NA. It is because dRTs are out of the range of forward dRT density distribution|
#| v. 2: Changes in how the function calculates dRTs. Avoiding NA PEPs. Marginally improves id rate. Only works with   |
#| libraries made by make.lib v.2. Cleaner output.                                                                     |
#+---------------------------------------------------------------------------------------------------------------------+

dRT <- function(path.in.lib, path.in.coeffs, path.in.evidence, path.out.PEP.updated){
  
  rt.lib <- read.delim(path.in.lib, header = TRUE, stringsAsFactors = FALSE)
  shift.coeffs <- read.delim(path.in.coeffs, header = TRUE, stringsAsFactors = FALSE)
  ev <- read.delim(path.in.evidence, header = TRUE, stringsAsFactors = FALSE)
  
  ev <- ev[, c('Raw.file', 'Sequence', 'PEP', 'Retention.time', "PIF", "Reporter.intensity.corrected.0",
               "Reporter.intensity.corrected.1", "Reporter.intensity.corrected.2", "Reporter.intensity.corrected.3",
               "Reporter.intensity.corrected.4", "Reporter.intensity.corrected.5", "Reporter.intensity.corrected.6",
               "Reporter.intensity.corrected.7", "Reporter.intensity.corrected.8", "Reporter.intensity.corrected.9",
               "id", "Leading.razor.protein")]
  ev <- ev[ev$Raw.file %in% shift.coeffs[!is.na(shift.coeffs$coeff), "experiment"], ]
  ev$protLabel <- substr(ev$Leading.razor.protein, 1, 3)
  
  # FORWARD:
  ev.for <- subset(ev, ev$protLabel=="sp|")
  ev.for$protLabel <- "FOR"
  ev.for <- ev.for[ev.for$Sequence %in% rt.lib$peptides, ]
  ev.for$RT.lib <- rt.lib$rt.median[match(ev.for$Sequence, rt.lib$peptides)]
  ev.for$intercept <- shift.coeffs$intercept[match(ev.for$Raw.file, shift.coeffs$experiment)]
  ev.for$coeff <- shift.coeffs$coeff[match(ev.for$Raw.file, shift.coeffs$experiment)]
  ev.for$RT.corrected <- ev.for$intercept + (ev.for$coeff * ev.for$Retention.time)
  ev.for$dRT.med <- abs(ev.for$RT.corrected-ev.for$RT.lib)
  for (i in 1:nrow(ev.for)){
    ev.for$dRT[i] <- min(abs(ev.for$RT.corrected[i]-na.omit(t(rt.lib[rt.lib$peptides==ev.for$Sequence[i], 2:(ncol(rt.lib)-4)]))))
  }
  ev.for <- ev.for[,-(20:21)]
  
  # REVERSE:
  ev.rev <- ev[ev$protLabel=="REV", ]
  set.seed(1)
  ev.rev$RT.lib <- sample(rt.lib$rt.median, nrow(ev.rev), replace=T)
  ev.rev$intercept <- shift.coeffs$intercept[match(ev.rev$Raw.file, shift.coeffs$experiment)]
  ev.rev$coeff <- shift.coeffs$coeff[match(ev.rev$Raw.file, shift.coeffs$experiment)]
  ev.rev$RT.corrected <- ev.rev$intercept + (ev.rev$coeff * ev.rev$Retention.time)
  ev.rev$dRT.med <- abs(ev.rev$RT.corrected - ev.rev$RT.lib)
  ev.rev$dRT <- ev.rev$dRT.med
  ev.rev <- ev.rev[,-(20:21)]
  
  ev.tot <- rbind(ev.rev, ev.for)
  row.names(ev.tot) <- NULL
  #ev.tot <- ev.perc.dRTed[sample(nrow(ev.perc.dRTed)), ]
  #row.names(ev.perc.dRTed) <- NULL
  
  # Updating PEPs:
  ev.PEP <- data.frame(Raw.file=character(), Sequence=character(), PEP=numeric(), Retention.time=numeric(), PIF=numeric(),
                       Reporter.intensity.corrected.0=numeric(), Reporter.intensity.corrected.1=numeric(),
                       Reporter.intensity.corrected.2=numeric(), Reporter.intensity.corrected.3=numeric(), 
                       Reporter.intensity.corrected.4=numeric(), Reporter.intensity.corrected.5=numeric(),
                       Reporter.intensity.corrected.6=numeric(), Reporter.intensity.corrected.7=numeric(),
                       Reporter.intensity.corrected.8=numeric(), Reporter.intensity.corrected.9=numeric(), id=numeric(),
                       Leading.razor.protein=character(), RT.lib=numeric(), RT.corrected=numeric(), dRT.med=numeric(),
                       dRT=numeric())
  exps <- unique(ev.tot$Raw.file)
  for (i in exps) {
    ex <- subset(ev.tot, ev.tot$Raw.file==i & ev.tot$PEP<=1)
    row.names(ex) <- NULL
    
    rev <- subset(ex, ex$protLabel=="REV")
    forw <- subset(ex, ex$protLabel=="FOR" & ex$PEP<0.02)
    
    den.for <- density(forw$dRT.med)
    den.rev <- density(rev$dRT.med)
    
    ex$Tr <- approx(den.for$x, den.for$y, xout=ex$dRT)$y
    ex$Fa <- approx(den.rev$x, den.rev$y, xout=ex$dRT)$y
    ex[is.na(ex$Tr), "Tr"] <- 1e-100
    ex$PEP.new <- (ex$PEP*ex$Fa)/((ex$PEP*ex$Fa)+((1-ex$PEP)*ex$Tr))
    
    ev.PEP <- rbind(ev.PEP, ex)
    
  }
  row.names(ev.PEP) <- NULL
  ev.PEP$Tr <- NULL
  ev.PEP$Fa <- NULL
  ev.PEP$dRT.med <- NULL
  
  write.table(ev.PEP, path.out.PEP.updated, sep = '\t', row.names = FALSE, quote = FALSE)
}








