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


