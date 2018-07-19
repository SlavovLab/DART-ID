library(tidyverse)
library(ggridges)
library(RColorBrewer)
source('Rscripts/lib.R')

## load data --------

source('Rscripts/alignment_comparison.R')

# or load data from file --------------------------------------------------
# generated 20180712
load('dat/error_df.RData')
load('dat/error_df_20180716.RData')

## or fake some data ------

df <- data.frame()
error_weights <- c(2, 3, 4, 2, 1, 2)
dat_len <- 300
#errors <- list()
for(i in 1:length(error_weights)) {
  pred <- runif(dat_len, min=5, max=60)
  obs <- pred + rnorm(dat_len, mean=0, sd=error_weights[i])
  df <- rbind(df, data.frame(
    pred=as.numeric(pred),
    obs=as.numeric(obs),
    error=as.numeric(obs-pred),
    type=as.character(i)
  ))
}

## scatterplots --------

get_density <- function(x, y, n = 100) {
  dens <- MASS::kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}
k <- 50
contour_cols <- viridis::viridis(k, alpha=0.3)
#cols <- brewer.pal(8,'Set1')[c(1:5, 8)]
cols <- c(brewer.pal(9, 'Blues')[5], brewer.pal(9, 'BuGn')[5], brewer.pal(9, 'Blues')[8],
          brewer.pal(9, 'Reds')[4], brewer.pal(9, 'PuRd')[5], brewer.pal(9, 'Reds')[7])

#pdf(file='manuscript/Figs/alignment_residual_scatter.pdf', width=4, height=6)
png(file='manuscript/Figs/alignment_residual_scatter_v5.png', width=3.5, height=5, units='in', res=250)

layout(rbind(c(1,4),c(3,5),c(2,6)))

par(oma=c(2.5, 2.5, 2, 0),
    mar=c(0.25, 0, 0.25, 0),
    cex.axis=1, pty='s')

for(i in 1:6) {
  
  #df_a <- df %>% filter(type==i)
  df_a <- error_df %>% filter(type == levels(error_df$type)[i])
  
  #dens <- get_density(df_a$pred_adjusted[!is.na(df_a$pred_adjusted)], df_a$obs[!is.na(df_a$pred_adjusted)], k)
  
  plot(df_a$pred_adjusted, df_a$obs, pch=16, cex=0.3, 
       #col=rgb(0,0,0,0.1),
       col=paste0(cols[i],'33'),
       #col=contour_cols[findInterval(dens, seq(0, max(dens), length.out=k))],
       xlab=NA, ylab=NA, xaxt='n', yaxt='n', xlim=c(5, 60), ylim=c(5, 60))
  abline(a=0, b=1, col='black')
  
  co <- cor(df_a$pred_adjusted, df_a$obs, use='pairwise.complete.obs')^2
  #sigma <- sd(df_a$error, na.rm=T)
  sigma <- mean(abs(df_a$error), na.rm=T)
  sigma <- formatC(sigma, digits=2, format='f')
  text(5, 56.5, 
       #paste0('p = ', formatC(co, digits=3)), 
       bquote(R^2*.(' = ')*.(formatC(co, digits=3, format='f'))),
       adj=c(0, 0.5), cex=1.5)
  #text(5, 52, 
  #     bquote(sigma*.(' = ')*.(sigma)*.(' min')),
  #     adj=c(0, 0.5), cex=1.1)
  #text(55, 10, paste(rep(i, 10), collapse=''), adj=c(1, 0.5))
  
  typ <- levels(error_df$type)[i]
  if(typ == 'MaxQuant MBR') typ <- 'MaxQuant'
  
  text(60, 7.5, typ, adj=c(1, 0.5), font=2, cex=1.5)
  
  if(i == 2 | i == 6) {
    axis(1, at=seq(10, 60, by=10), tck=-0.02, mgp=c(0, 0.25, 0))
  } else {
    axis(1, at=seq(10, 60, by=10), tck=-0.02, labels=NA)
  }
  if(i <= 3) {
    axis(2, at=seq(10, 60, by=10), tck=-0.02, mgp=c(0, 0.4, 0), las=1)
  } else {
    axis(2, at=seq(10, 60, by=10), tck=-0.02, labels=NA)
  }
}

mtext('Observed RT (min)', outer=T, side=2, line=1.2, cex=1)
mtext('Predicted RT (min)', outer=T, side=1, line=1.4, cex=0.9, adj=0.1)
mtext('Aligned RT (min)', outer=T, side=1, line=1.4, cex=0.9, adj=0.9)
mtext('Prediction', outer=T, side=3, line=0.1, cex=1, adj=0.18, font=2)
mtext('Alignment', outer=T, side=3, line=0.1, cex=1, adj=0.83, font=2)

dev.off()


## ridgeplots -------

cols <- c(brewer.pal(9, 'Blues')[5], brewer.pal(9, 'BuGn')[5], brewer.pal(9, 'Blues')[8],
          brewer.pal(9, 'Reds')[4], brewer.pal(9, 'PuRd')[5], brewer.pal(9, 'Reds')[7])

df_a <- error_df %>% filter(type %in% c('SSRCalc', 'ELUDE', 'BioLCCC'))
errors <- df_a %>% group_by(type) %>% summarise(err=mean(abs(error), na.rm=T)) %>% pull(err)
p <- 
ggplot(df_a) +
  geom_density_ridges(aes(error, type, fill=type),
                      #stat='binline', bins=51, 
                      size=0.5, color=NA,
                      rel_min_height=0.01, scale=1.25) +
  scale_x_continuous(limits=c(-10, 10), breaks=seq(-10, 10, by=5)) +
  scale_y_discrete(limits=c('ELUDE', 'BioLCCC', 'SSRCalc'), 
                   #labels=paste0(c('ELUDE\n', 'SSRCalc\n', 'BioLCCC\n'), 
                    #             #parse(text=c('sigma', 'sigma', 'sigma')),
                    #             c(bquote(sigma), bquote(sigma), bquote(sigma)),
                    #            formatC(errors, digits=3, format='f')),
                   #labels=c(bquote(atop(ELUDE, sigma*.('=')*.(formatC(errors[2], digits=2, format='f'))*.(' min'))),
                  #          bquote(atop(BioLCCC, sigma*.('=')*.(formatC(errors[3], digits=2, format='f'))*.(' min'))),
                   #         bquote(atop(SSRCalc, sigma*.('=')*.(formatC(errors[1], digits=2, format='f'))*.(' min')))
                  #          ),
                  labels=c('ELUDE', 'BioLCCC', 'SSRCalc'),
                   expand=c(0.02,0)) +
  scale_fill_manual(values=cols, guide=F) +
  labs(x='Residual RT (min)', y=NULL) +
  theme_ridges() + #theme_bert()
  theme(
    axis.title.x=element_text(hjust=0.5, size=10),
    axis.text=element_text(size=8),
    axis.title.y=element_text(hjust=0.5)
  )

ggsave('manuscript/Figs/alignment_binline_pred_v2.pdf', p, device='pdf', width=1.75, height=2.5, units='in')

df_b <- error_df %>% filter(type %in% c('iRT', 'MaxQuant MBR', 'DART-ID'))
errors <- df_b %>% group_by(type) %>% summarise(err=mean(abs(error), na.rm=T)) %>% pull(err)
p <- 
  ggplot(df_b) +
  geom_density_ridges(aes(error, type, fill=type),
                      #stat='binline', bins=51, 
                      size=0.5, color=NA, panel_scaling=F,
                      rel_min_height=0.01, scale=1.25) +
  scale_x_continuous(limits=c(-2, 2), breaks=seq(-2, 2, by=1)) +
  scale_y_discrete(limits=c('DART-ID', 'MaxQuant MBR', 'iRT'), expand=c(0.02,0),
                   #labels=c(bquote(atop(DART*.('-')*ID, sigma*.('=')*.(formatC(errors[3], digits=2, format='f'))*.(' min'))),
                    #        bquote(atop(MaxQuant, sigma*.('=')*.(formatC(errors[2], digits=2, format='f'))*.(' min'))),
                     #       bquote(atop(iRT, sigma*.('=')*.(formatC(errors[1], digits=2, format='f'))*.(' min')))
                   #)
                   labels=c('DART-ID', 'MaxQuant', 'iRT')
                   ) +
  scale_fill_manual(values=cols[4:6], guide=F) +
  labs(x='Residual RT (min)', y=NULL) +
  theme_ridges() + #theme_bert()
  theme(
    plot.margin=margin(1, 0.2, 0.2, 0.2, 'cm'),
    axis.title.x=element_text(hjust=0.5, size=10),
    axis.text=element_text(size=8),
    axis.title.y=element_text(hjust=0.5)
  )

ggsave('manuscript/Figs/alignment_binline_align_v2.pdf', p, device='pdf', width=1.75, height=2.5, units='in')

# boxplots -----------------------------------------------------------------

types <- levels(error_df$type)

png(file='manuscript/Figs/alignment_residual_boxplot_v3.png', width=3.5, height=2.5, units='in', res=250)

par(mar=c(2.75, 2.5, 0.25, 0.25))

#boxplot(error~type, data=error_df, horizontal=T)
boxplot(error~type, data=error_df, #horizontal=T,
        xlab=NA, ylab=NA, xaxt='n', yaxt='n',
        col=cols,
        ylim=c(-10, 10),
        range=1.1,
        outpch=16, outcex=0.5, outcol=rgb(0,0,0,0.1))

axis(2, 
     #at=seq(-10, 10, by=4), 
     at=c(-10, -6, -3, 0, 3, 6, 10),
     tck=-0.02,
     mgp=c(0, 0.5, 0), las=1)
#axis(1, at=seq(1,6), labels=types, tck=-0.02,
#     mgp=c(0, 0.5, 0), las=2)
axis(1, at=seq(1,6), labels=NA, tck=-0.02)
text(seq(1, 6), par('usr')[3]-1.1, srt=25, adj=c(1, 0.5), xpd=T,
     labels=types, cex=0.75)

mtext('Residual RT (min)', side=2, line=1.5, cex=1)

dev.off()


