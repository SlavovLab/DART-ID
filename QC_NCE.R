library(tidyverse)
library(pracma)
library(reshape2)
source('lib.R')

ev <- read_tsv("~/git/RTLib/Alignments/NCE_20180520_5/ev_combined.txt")

exp <- gsub('[0-9]{6}[A-Z]{1}_QC_SQC', '', ev$`Raw file`)
# should be 53 levels here
exp_factor <- as.numeric(as.factor(exp))

collision_energies <- c(30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, #SQC55A-K
                        #3040, 3045, 3050, #SQC55L-N - stepped
                        40, 45, 50, #SQC55L-N - stepped
                        38,39,40,41,42,43,44,45,46,47, #SQC57A-J
                        44,44, #SQC61A-B
                        48,47,46,45,44,43,42,41,40,39,38,37,36, #SQC62A-M
                        37,37,37, #SQC65A1-3
                        47,47,47, #SQC65B1-3
                        37,37, #SQC65C1-2
                        47,47, #SQC65D1-2
                        37,37, #SQC65E1-2
                        47,47) #SQC65F1-2

ev$collision_energies <- collision_energies[exp_factor]
ev$dPEP <- -log10(ev$pep_new / ev$PEP)

ev_f <- ev %>% sample_n(1e4)
plot(ev_f$PEP, ev_f$pep_new, log='xy', xlim=c(1e-10, 1), ylim=c(1e-10, 1), pch=16, cex=0.2)
abline(a=0, b=1, col='red')

plot(ev_f$collision_energies, ev_f$dPEP, pch=16, cex=0.2)

ggplot(ev) +
  geom_violin(aes(x=as.factor(collision_energies), y=dPEP), scale='area', fill='black')

pep_threshold <- 0.01

id_rates <- ev %>%
  filter(!grepl('SQC55[L-N]{1}', `Raw file`)) %>%
  mutate(collision_energies=as.factor(collision_energies)) %>%
  
  group_by(collision_energies) %>%
  summarise(identified=sum(PEP < pep_threshold) / length(unique(`Raw file`)), 
            newly_identified=sum(pep_updated < pep_threshold, na.rm=T) / length(unique(`Raw file`)),
            all=n() / length(unique(`Raw file`)))

b <- ggplot(melt(id_rates)) +
  geom_bar(aes(x=collision_energies, y=value, fill=variable), stat='identity', position='dodge') +
  scale_y_continuous(expand=c(0,0)) +
  labs(x="Collision Energy", y="Average PSMs per Run", fill="",
       title=paste("ID Efficiency by Collision Energy\nPEP Threshold:", pep_threshold)) +
  theme_bert()

# relative increase
id_rates %>%
  mutate(rel_increase=newly_identified/identified) %>%
  ggplot() +
  geom_bar(aes(x=collision_energies, y=rel_increase), stat='identity') +
  labs(x="Collision Energy", y=paste("Relative Increase of PSMs\nIDed with PEP <", pep_threshold))
  
# stepped vs. unstepped
ev$stepped <- "Unstepped"
ev$stepped[grepl('SQC55[L-N]{1}', ev$`Raw file`)] <- "Stepped"

id_rates <- ev %>%
  group_by(stepped) %>%
  summarise(identified=sum(PEP < 0.05) / length(unique(`Raw file`)), 
            newly_identified=sum(pep_updated < 0.05, na.rm=T) / length(unique(`Raw file`)),
            all=length(PEP) / length(unique(`Raw file`)))

ggplot(melt(id_rates)) +
  geom_bar(aes(x=stepped, y=value, fill=variable), stat='identity', position='dodge')

# scatter
ev %>% sample_n(1e4) %>%
  ggplot() +
  geom_point(aes(x=collision_energies, y=dPEP))

# look at RI Intensities
ev %>%
  mutate(RI.mean=`Reporter intensity corrected 0` + `Reporter intensity corrected 1` + `Reporter intensity corrected 2` + `Reporter intensity corrected 3` + `Reporter intensity corrected 4` + `Reporter intensity corrected 5` + `Reporter intensity corrected 6` + `Reporter intensity corrected 7` + `Reporter intensity corrected 8` + `Reporter intensity corrected 9`) %>%
  mutate(RI.mean=sapply(RI.mean, function(x) { if(x==0) return(NA) else return(x) })) %>%
  group_by(collision_energies) %>%
  summarise(ri.intensity=mean(RI.mean, na.rm=T)) %>%
  ggplot() +
  geom_bar(aes(x=collision_energies, y=ri.intensity), stat='identity') +
  scale_x_continuous(breaks=seq(30, 50, 1), labels=seq(30, 50, 1)) +
  labs(x="Collision Energy", y="mean(RI Intensity)", fill="Channel")

# individual RI intensities
ev %>%
  group_by(collision_energies) %>%
  summarise(
    `126C`=mean(`Reporter intensity corrected 0`),
    `127N`=mean(`Reporter intensity corrected 1`),
    `127C`=mean(`Reporter intensity corrected 2`),
    `128N`=mean(`Reporter intensity corrected 3`),
    `128C`=mean(`Reporter intensity corrected 4`),
    `129N`=mean(`Reporter intensity corrected 5`),
    `129C`=mean(`Reporter intensity corrected 6`),
    `130N`=mean(`Reporter intensity corrected 7`),
    `130C`=mean(`Reporter intensity corrected 8`),
    `131N`=mean(`Reporter intensity corrected 9`)
  ) %>%
  gather(variable, value, -collision_energies) %>%
  ggplot() +
  geom_bar(aes(x=collision_energies, y=log(value), fill=variable), stat='identity') +
  facet_grid(variable~.) +
  scale_x_continuous(breaks=seq(30, 50, by=1)) +
  scale_y_continuous(breaks=c(0, 5, 10)) +
  labs(x="Collision Energy", y="log(mean(RI Intensity))", fill="Channel")
