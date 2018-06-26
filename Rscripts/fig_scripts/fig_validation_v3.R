library(tidyverse)
source('Rscripts/lib.R')
source('Rscripts/validate.lib.3.R')

if(!exists("cvs_all") | is.null(cvs_all)) {
  cvs_all <- validate.lib.3(ev)
}

# ECDF?
# fig2d <-
# ggplot(cvs_all) +
#   stat_ecdf(aes(value, color=Method), geom='step', position='identity', size=0.75) +
#   scale_color_manual(values=c(av[1], av[2], av[3], av[5])) +
#   scale_x_continuous(limits=c(0, 0.4)) +
#   labs(x="CV of Quantitation", y="Cumulative Density", color=NULL,
#        title="Quantitation consistency of\nPSMs within proteins") +
#   guides(color=guide_legend(nrow=4, override.aes=list(size=1.5))) +
#   theme_bert() + theme(
#     plot.margin=margin(0.1,0.3,0.1,0.1,'cm'),
#     axis.text.y=element_blank(),
#     legend.position=c(0.33, 0.88),
#     legend.key.height=unit(0.3, 'cm'),
#     legend.box.spacing=unit(c(0,0,0,0),'cm'),
#     plot.title=element_text(size=9, lineheight=0.9, margin=margin(0,0,0.1,0,'cm'))
#   )

cvs_all$Method <- factor(cvs_all$Method, labels=c("Spectra", "DART-ID", "Percolator", "Null"))
#cvs_all$Method <- factor(cvs_all$Method, labels=c("Spec", "Spec+RT", "Perc", "Null"))

fig2d <-
  ggplot(cvs_all) +
  geom_boxplot(aes(Method,value, group=Method, fill=Method),color='black',
               outlier.shape='x', outlier.size=3, outlier.alpha=0) +
  #geom_violin(aes(Method, value, group=Method, fill=Method)) +
  scale_y_continuous(limits=c(0.175, 0.35)) +
  scale_fill_manual(values=c(av[1], av[2], av[3], av[5]), guide=F) +
  labs(x=NULL, y="CV of Relative Quantitation", color=NULL,
       title="Consistency of\nProtein Quantitation") +
  theme_bert() + theme(
    plot.margin=margin(0.1,0.3,0.1,0.1,'cm'),
    plot.title=element_text(size=10, margin=margin(0,0,0.1,0,'cm'), hjust=0.5),
    axis.title=element_text(size=10),
    axis.text.x=element_text(size=10, angle=30, hjust=1, vjust=1)
    #axis.text.x=element_text(size=8)
  )

# t-tests
# t.test(cvs_all$value[cvs_all$Method=="Spectra"], 
#        cvs_all$value[cvs_all$Method=="Null"], var.equal=T)
# t.test(cvs_all$value[cvs_all$Method=="DART-ID"], cvs_all$value[cvs_all$Method=="Null"], var.equal=T)
# t.test(cvs_all$value[cvs_all$Method=="Percolator"], cvs_all$value[cvs_all$Method=="Null"], var.equal=T)
# 
# var.test(cvs_all$value[cvs_all$Method=="Spectra"], 
#          cvs_all$value[cvs_all$Method=="DART-ID"])
# t.test(cvs_all$value[cvs_all$Method=="Spectra"], 
#        cvs_all$value[cvs_all$Method=="DART-ID"], var.equal=T)


#return(fig2d)
