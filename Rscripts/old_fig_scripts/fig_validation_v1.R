## heatmap ------

save(prot.data, prot.names, prot.cor, prot.data.new, prot.names.new, prot.cor.new, file='Fig2dat.RData')

# this is data pulled from the guts of the validate.lib() function.
# options run were:
# ev = dat/ev.adj.Fit3c.txt
# raw.files = c("160822A_NC_set#19A_180min_200NL_60C+N=3_250ms", "160822A_NC_set#19A_180min_200NL_60C+N=3_750ms", "160822A_NC_set#19A_180min_200NL_60C+N=3_500ms")
# remove.REV = T
# remove.CON = T
# remove.keratins = T
# sparse.filter = 0.5
# unique.only = F
# norm.cols = T
# norm.rows = T
# psm.freq.thresh = 10
# alpha = 0.05
# cor.size.thresh = 5

# even though we already have this nice correlation matrix
# we need to get it in long form for ggplot

df <- melt(prot.cor.new, varnames=c('Col', 'Row'), value.name='Cor')
df %>%
ggplot(aes(x=Col, y=Row, fill=Cor)) +
  geom_tile() +
  scale_x_continuous(limits=c(0+0.5, max(df$Col)+0.5),
                     breaks=seq(1,max(df$Col)),
                     labels=prot.names.new,
                     expand=c(0,0)) +
  scale_y_continuous(limits=c(0+0.5, max(df$Row)+0.5),
                     breaks=seq(1,max(df$Row)),
                     labels=prot.names.new,
                     expand=c(0,0.1)) +
  #scale_fill_gradient2(low='red', mid='black', high='green', 
  #                     midpoint=0, na.value='white') +
  scale_fill_gradientn(colors=brewer.pal('RdBu', n=11), na.value=brewer.pal('RdBu', n=11)[11]) +
  labs(x=NULL, y=NULL) +
  theme_bw(base_size=16, base_family="Helvetica") %+replace% 
  theme(
    panel.background=element_blank(),
    panel.border=element_blank(),
    #panel.spacing=unit(0, 'in'),
    panel.grid=element_blank(),
    axis.text=element_text(size=10),
    axis.text.x=element_text(angle=45, hjust=1, vjust=1),
    axis.ticks=element_blank()
    #plot.margin=margin(t=0,r=0,b=0,l=0)
  )


## group cor density plot -----


## Group plots -----

library(gridExtra)
library(grid)
library(gtable)

exps <- c('160712A_NC_19A[0-9]{1}', '160801A_NC_set19A',
          '160818A_NC_set#19A', '160822A_NC_set#19A')
g <- gtable()

cor.plots <- list()
for(i in 1:length(exps)) {
  
  # TODO:
  # color with a palette from RColorBrewer
  # get cowplot and use that?
  # dont plot here, get density and store that in a df with the experiment name. 
  # we'll use the default faceting later
  
  raw.files <- all.raw.files[grep(exps[i], all.raw.files)]
  cors <- validate.lib(ev, raw.files=raw.files)
  
  p <- ggplot(cors, aes(x=data, color=Type)) +
  #ggplot(cors, aes(x=data, color=Type)) +
    geom_line(stat='density', position='identity') +
    scale_color_manual(values=c(
      "New"="blue",
      "Original"="black",
      "Null"="red"
    ), guide=F) +
    scale_x_continuous(breaks=c(-1, -0.5, 0, 0.5, 1)) +
    scale_y_continuous(limits=c(0, 2), breaks=c(0, 0.5, 1, 1.5, 2)) +
    labs(title=NULL, x=NULL, y=NULL, color=NULL, fill=NULL) +
    theme_bert()
  
  
  
  pGrob <- ggplotGrob(p)
  
  cor.plots <- c(cor.plots, list(pGrob))
}

grid.arrange(cor.plots[[1]], cor.plots[[2]],
             cor.plots[[3]], cor.plots[[4]],
             nrow=2, ncol=2)


p <- ggplot(cors, aes(x=data, color=Type)) +
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

ggsave(filename=paste0('validation_figs/',i,'.pdf'), plot=p, device='pdf',
       width=5, height=5, units='in')