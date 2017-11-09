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

dRT <- function(path.in.lib='dat/RT.lib.elite.txt', 
                path.in.coeffs='dat/shift.coeffs.elite.txt', 
                path.in.evidence='dat/evidence.txt', 
                path.out.PEP.updated='dat/evidence+dRT.elite.txt'){
  
  rt.lib <- read.delim(path.in.lib, header = TRUE, stringsAsFactors = FALSE)
  shift.coeffs <- read.delim(path.in.coeffs, header = TRUE, stringsAsFactors = FALSE)
  ev <- read.delim(path.in.evidence, header = TRUE, stringsAsFactors = FALSE)
  
  ev <- ev[, c('Raw.file', 'Sequence', 'PEP', 'Retention.time', "PIF", 
               "Reporter.intensity.corrected.0", "Reporter.intensity.corrected.1", 
               "Reporter.intensity.corrected.2", "Reporter.intensity.corrected.3",
               "Reporter.intensity.corrected.4", "Reporter.intensity.corrected.5", 
               "Reporter.intensity.corrected.6", "Reporter.intensity.corrected.7", 
               "Reporter.intensity.corrected.8", "Reporter.intensity.corrected.9",
               "id", "Leading.razor.protein", "Proteins")]
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
  #for (i in 1:nrow(ev.for)){
    # smallest difference between corrected RT and any one of the RT entries in the
    # RT library
  #  ev.for$dRT[i] <- min(
  #    abs(ev.for$RT.corrected[i] -
  #          na.omit(t(rt.lib[rt.lib$peptides==ev.for$Sequence[i], 
  #                   2:(ncol(rt.lib)-4)]))
  #    )
  #  )
    # this takes a while
  #  flush.console()
  #  cat('\r', i, '/', nrow(ev.for), '              ')
  #}
  ev.for <- ev.for[,-c(21,22)]
  
  # REVERSE:
  ev.rev <- ev[ev$protLabel=="REV", ]
  set.seed(1)
  ev.rev$RT.lib <- sample(rt.lib$rt.median, nrow(ev.rev), replace=T)
  ev.rev$intercept <- shift.coeffs$intercept[match(ev.rev$Raw.file, shift.coeffs$experiment)]
  ev.rev$coeff <- shift.coeffs$coeff[match(ev.rev$Raw.file, shift.coeffs$experiment)]
  ev.rev$RT.corrected <- ev.rev$intercept + (ev.rev$coeff * ev.rev$Retention.time)
  ev.rev$dRT.med <- abs(ev.rev$RT.corrected - ev.rev$RT.lib)
  #ev.rev$dRT <- ev.rev$dRT.med
  ev.rev <- ev.rev[,-c(21,22)]
  
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
                       Leading.razor.protein=character(), Proteins=character(), RT.lib=numeric(), RT.corrected=numeric(), dRT.med=numeric(),
                       dRT=numeric())
  
  dens.forw <- list()
  dens.rev <- list()
  dRT.forw <- list()
  dRT.rev <- list()
  
  exps <- unique(ev.tot$Raw.file)
  counter <- 0
  
  for (i in exps) {
    counter <- counter + 1
    
    cat('\r', 'Processing ', counter, '/', length(exps), ' ', i,  '                                          ')
    flush.console()
    
    ex <- subset(ev.tot, ev.tot$Raw.file==i & ev.tot$PEP<=1)
    row.names(ex) <- NULL
    
    rev <- subset(ex, ex$protLabel=="REV")
    forw <- subset(ex, ex$protLabel=="FOR" & ex$PEP<0.02)
    
    den.for <- density(forw$dRT.med)
    den.rev <- density(rev$dRT.med)
    
    dens.forw[[counter]] <- den.for
    dens.rev[[counter]] <- den.rev
    #dRT.forw[[counter]] <- forw$dRT.med
    #dRT.rev[[counter]] <- rev$dRT.med
    
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

#par(mfrow=c(2,2))
#for(i in sample.int(length(exps), size=4)) {
#  plot(dens.forw[[i]], col='blue', xlab='dRT (min)', 
#       main=paste0('Correct (blue) vs. Incorrect (red)\n', 'Exp: ', i))
#  lines(dens.rev[[i]], col='red')
#}

