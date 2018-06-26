validate.lib.3 <- function(ev, exclude.cols=c(1:4, 7:11)) {
  library(reshape2)
  
  # find out the indices of columns with reporter ion data
  # should be 10 columns, but the code doesn't rely on this number
  data.cols <- grep('Reporter intensity corrected', colnames(ev))
  cat("Parsed columns '", paste(colnames(ev)[data.cols], collapse=", "), "' as RI data columns\n")
  
  # ignore empty, carrier channels
  cat("Excluding columns", exclude.cols, "\n")
  data.cols <- data.cols[-exclude.cols]
  cat("Remaining columns: '", paste(colnames(ev)[data.cols], collapse=", "), "' as RI data columns\n")

  cat("Processing proteins, and removing contaminants and decoys...\n")  
  # extract uniprot ID
  ev.f <- ev %>%
    mutate(Protein=sapply(strsplit(`Leading razor protein`, "\\|"), function(p) {
      if(length(unlist(p)) == 1) return(p[1])
      else if(length(unlist(p)) == 3) return(p[2])
      else return(p[1])
    })) %>%
    filter(!grepl("REV__", Protein)) %>%
    filter(!grepl("CON__", Proteins)) %>%
    # only take SQC experiments
    filter(grepl("SQC", `Raw file`)) %>%
    filter(!is.na(Protein))
  
  cat("Removing proteins without single-cell RI data\n")
  ev.f <- ev.f %>% filter(apply(ev.f[,data.cols]!=0, 1, sum, na.rm=T) == length(data.cols))
  
  cat(nrow(ev.f), "observations remaining\n")
  
  ev.f[,data.cols] <- data.matrix(ev.f %>% select(data.cols))
  
  # set zeroes to NA
  ev.f[,data.cols][ev.f[,data.cols]==0] <- NA
  
  ## normalize data
  
  cat("Normalizing data...\n")
  # first normalize by column, by median
  ev.f[,data.cols] <- t(t(ev.f[,data.cols]) / apply(ev.f[,data.cols], 2, median, na.rm=T))
  # now normalize across rows, by mean
  ev.f[,data.cols] <- ev.f[,data.cols] / apply(ev.f[,data.cols], 1, mean, na.rm=T)
  
  # this data has 6 columns, with alternating J, U cells
  # select peptides with some solid signal, and differences between the two cell types
  
  
  pep_thresh <- 0.01
  prot.psm.thresh <- 50
  
  prots_orig <- ev.f %>% 
    filter(PEP < pep_thresh) %>% 
    group_by(Protein) %>%
    summarise(n=n()) %>%
    filter(n > prot.psm.thresh) %>%
    arrange(desc(n)) %>%
    pull(Protein)
  
  prots_new <- ev.f %>%
    filter(PEP > pep_thresh & pep_new < pep_thresh) %>%
    group_by(Protein) %>%
    summarise(n=n()) %>%
    filter(n > prot.psm.thresh) %>%
    arrange(desc(n)) %>%
    pull(Protein)
  
  # only take the intersection with original proteins
  prots_new <- prots_new[prots_new %in% prots_orig]
  
  prots <- prots_new
  
  prot_cvs_orig <- zeros(length(prots), length(data.cols))
  prot_cvs_new  <- zeros(length(prots), length(data.cols))
  prot_cvs_perc <- zeros(length(prots), length(data.cols))
  prot_cvs_null <- zeros(length(prots), length(data.cols))
  
  for(i in 1:length(prots)) {
    cat('\r', i, '/', length(prots), '-', prots[i], '                           ')
    flush.console()
    
    ev.a <- ev.f %>% filter(Protein==prots[i])
    dmat <- data.matrix(ev.a %>% filter(PEP < pep_thresh) %>% select(data.cols))
    prot_cvs_orig[i,] <- apply(dmat, 2, sd) / apply(dmat, 2, mean)
    dmat <- data.matrix(ev.a %>% filter(PEP > pep_thresh & pep_new < pep_thresh) %>% select(data.cols))
    prot_cvs_new[i,] <- apply(dmat, 2, sd) / apply(dmat, 2, mean)
    dmat <- data.matrix(ev.a %>% filter(PEP > pep_thresh & pep_perc < pep_thresh) %>% select(data.cols))
    prot_cvs_perc[i,] <- apply(dmat, 2, sd) / apply(dmat, 2, mean)
    
    # generate null distribution
    dmat <- data.matrix(ev.f %>% filter(PEP < pep_thresh) %>% sample_n(50) %>% select(data.cols))
    prot_cvs_null[i,] <- apply(dmat, 2, sd) / apply(dmat, 2, mean)
  }
  
  cat("\nCollecing results...\n")
  cvs_all <- data.frame()
  cvs_all <- rbind(cvs_all, melt(prot_cvs_orig)  %>% mutate(Method="Spectra"))
  cvs_all <- rbind(cvs_all, melt(prot_cvs_new)   %>% mutate(Method="Spectra+RT"))
  cvs_all <- rbind(cvs_all, melt(prot_cvs_perc)  %>% mutate(Method="Percolator"))
  cvs_all <- rbind(cvs_all, melt(prot_cvs_null)  %>% mutate(Method="Null"))
  
  cvs_all$Method <- factor(cvs_all$Method, levels=c("Spectra", "Spectra+RT", "Percolator", "Null"))
  #factor(cvs_all$Var2, labels=c("", "b", "c", "d", "e", "f"))
  
  return(cvs_all)
}
