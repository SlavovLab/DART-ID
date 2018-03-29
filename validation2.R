library(tidyverse)
source('lib.R')

ev <- read_tsv('dat/ev.adj.Fit3c.RTL.txt')

all.raw.files <- unique(ev$`Raw file`)
files.19A <- all.raw.files[grep('19A', all.raw.files)]

# lets take only the files run on 160712 (7/16/12)
raw.files <- all.raw.files[grep('160712A_NC_19A', all.raw.files)]

raw.files <- all.raw.files[grep('38A', all.raw.files)]

#exps <- c('160712A_NC_19A','38A')
#exps <- c('38B', '38C', '38D')
#exps <- c('18A1', '18A2', '18B1', '18B2', '18C1', '18C2')
#exps <- all.raw.files

exps <- c('set#37A', 'set#37B', 'set#37C', 'set#37D',
          'set#38A', 'set#38B', 'set#38C', 'set#38D', 'set#39A', 'set#39B',
          'set#40A', 'set#40B', 'set#40C', 'set#40D',
          'set#32B', 'set#32C', 'set#32D',
          paste0('set#30', LETTERS[1:10]), 'set#29B',
          'set#29A', 'set#29C',
          paste0('set24', c('A','B','C')),
          paste0('set25', c('A','B','C')),
          paste0('set#26', c('A','B','C','D','E')),
          paste0('set#27', c('A','B','C','D')),
          paste0('set#28', c('A','B','C','D')),
          'set19B', '160712A_NC_19A',
          'set#30', 'set#29', '160822A_NC_set#27', '160913A_NC_set#32',
          '160712A_NC_19A', '160720A_NC_19A', '160727A_NC_19A', 
          '160801A_NC_set19A', '160803A_NC_set19A', '160808A_NC_set19A',
          '160818A_NC_set#19A', '160822A_NC_set#19A', '160831A_NC_set#19A',
          '160927A_NC_set#19A')
exps <- c('160712A_NC_19A')

for(i in exps) {
  #pdf(paste0('validation_figs/',i,'.pdf'), width=5, height=5)
  #cat(i,'\n')
  raw.files <- all.raw.files[grep(i, all.raw.files)]
  #raw.files <- i
  
  cors <- NULL
  tryCatch({
    #cors <- validate.lib(ev, raw.files=raw.files, unique.only=T, psm.freq.thresh=20)
    cors <- validate.lib(ev, raw.files=raw.files)
  }, error=function(e) {
    print(e)
  })
  if(is.null(cors)) next
  
  #cors <- cors %>%
  #  filter(Type != 'Total')
  
  
  
  #p <- 
  ggplot(cors, aes(x=data, color=Type)) +
    #geom_histogram(aes(y=..density.., color=NULL, fill=Type), position='identity', alpha=0.5, bins=30) +
    geom_line(stat='density', position='identity') +
    scale_color_manual(values=c(
      "New"="blue",
      "Original"="black",
      "Null"="red"
    )) +
    scale_fill_manual(values=c(
      "New"="blue",
      "Original"="black",
      "Null"="red"
    ), guide=FALSE) +
    labs(title=paste(i, 'Validation'),
         x='Correlations', y='Density', color=NULL, fill=NULL) +
    theme_bert()
  
  #ggsave(filename=paste0('validation_figs/',i,'.pdf'), plot=p, device='pdf',
  #       width=5, height=5, units='in')
  
  #dev.off()
}

data <- data.matrix(ev[,data.cols])
# look at distributions of data in different RI channels
data <- t(t(data) / apply(data, 2, median, na.rm=T))
data <- data / apply(data, 1, median, na.rm=T)
#data <- data / apply(data, 1, mean, na.rm=T)

# melt into df
data <- melt(data)
colnames(data) <- c('PSM', 'Channel', 'Intensity')
# remove NAs
data <- data[!is.na(data$Intensity),]

levels(data$Channel) <- c('126', '127N', '127C', '128N', '128C', '129N', '129C', '130N', '130C', '131')

ggplot(data, aes(x=Intensity, y=..density.., fill=Channel)) +
  geom_density() +
  #geom_histogram(bins=50) +
  facet_grid(Channel~.) +
  scale_x_continuous(limits=c(0, 3)) +
  #scale_x_continuous(limits=c(0, 2000)) +
  labs(title='Normalized RI Levels - Row Median')

ggplot(data) +
  geom_boxplot(aes(x=Channel, y=Intensity), outlier.size = 0) +
  scale_y_continuous(limits=c(0,10))

colors <- brewer.pal('Set3', n=10)
df <- data.frame()
for(i in 1:10) {
  b <- density(data[,i], na.rm=T)
  df <- rbind(df, data.frame(
    x=as.numeric(b$x),
    y=as.numeric(b$y),
    channel=as.character(rep(i, length(b$x)))
  ))
}

ggplot(df, aes(x=x, y=y, color=channel)) +
  geom_path() +
  scale_x_continuous(limits=c(0,10)) +
  scale_color_manual(values=colors)
