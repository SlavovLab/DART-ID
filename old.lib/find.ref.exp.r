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
  counter <- 0
  for (ii in exps.all) {
    C10 <- 0 
    counter <- counter + 1
    
    cat('\r', 'Processing ', counter, '/', length(exps.all), ' ', i,  '                                          ')
    flush.console()
    
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



