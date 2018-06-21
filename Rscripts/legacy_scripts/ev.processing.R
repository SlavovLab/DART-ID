ev <- read_tsv('dat/ev.adj.elite.txt')

ev.f <- ev[,c("Sequence", "Raw file", "Leading razor protein",
                "Reporter intensity corrected 0", "Reporter intensity corrected 1", 
                "Reporter intensity corrected 2", "Reporter intensity corrected 3",
                "Reporter intensity corrected 4", "Reporter intensity corrected 5", 
                "Reporter intensity corrected 6", "Reporter intensity corrected 7", 
                "Reporter intensity corrected 8", "Reporter intensity corrected 9")]
colnames(ev.f)[c(2,3)] <- c('Raw.file', 'Protein')
ev.f$ID <- seq(1,nrow(ev.f))

ev.dat <- as.matrix(ev.f[,-c(1,2,3)])
colnames(ev.dat) <- NULL
# 0 = NA
ev.dat[ev.dat == 0] <- NA
# invalidate empty and carrier channel
ev.dat[,c(8, 10)] <- NA

# remove rows that have no quantitation
no.quant <- apply(ev.dat, MARGIN=1, FUN=sum, na.rm=TRUE) == 0
ev.f <- ev.f[!no.quant,]

# melt ev.f
ev.m <- melt(ev.f, id.vars=c('Sequence', 'Raw.file', 'Protein', 'ID'), na.rm=FALSE,
             variable.name='Channel', value.name='Intensity')
# more sensible channel names
ev.m$Channel <- as.numeric(ev.m$Channel)
# exclude empty and carrier channel
ev.m <- subset(ev.m, !(Channel %in% c(8, 10)))
# 0 = NA
ev.m$Intensity[ev.m$Intensity == 0] <- NA

# filtering
# remove CON and REV proteins
prot.label <- substr(ev.m$Protein, 1, 3)
ev.m <- ev.m[-grep('CON|REV', prot.label),]

# get peptide levels - collapse PSMs from same peptide, raw file, and channel
ev.pep <- aggregate(Intensity ~ Sequence + Raw.file + Channel, data=ev.m, FUN=median)
# get protein levels - collapse PSMs from same raw file, channel
ev.prot <- aggregate(Intensity ~ Protein + Raw.file + Channel, data=ev.m, FUN=median)

# protein count
ev.prot.freq <- aggregate(.~Protein+Raw.file+Channel, data=ev.m, FUN=sum)
 