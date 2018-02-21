filter.exps <- function(ev, pep.thresh=1e-2) {
  # ev <- read_tsv('dat/evidence_elite.txt')
  
  #removed.exps <- c()
  
  # confidence threshold - harsh to increase confidence
  # and decrease computation time
  # pep.thresh <- 1e-2
  
  # make sure we have enough observations in all experiments
  # lets say... need at least 100 to get a good idea of experiment quality
  #below.thresh <- ev %>% 
  #  group_by(`Raw file`) %>% 
  #  summarise(n=sum(PEP < pep.thresh)) %>%
  #  arrange(desc(n)) %>%
  #  filter(n < 100) %>%
  #  pull(`Raw file`)
  
  #ev <- ev %>% filter(!`Raw file` %in% below.thresh)
  
  raw.files <- unique(ev$`Raw file`)
  cor.mat <- matrix(nrow=length(raw.files), ncol=length(raw.files))
  
  # correlate high-confidence PSMs between experiments
  # only need the upper triangular of the correlation matrix
  for(i in 1:length(raw.files)) {
    ev.a <- ev %>% 
      filter(PEP < pep.thresh) %>%
      filter(`Raw file`==raw.files[i]) %>%
      select(c('Mod. peptide ID', 'Retention time'))
    
    if(i == 1:length(raw.files)) next
    
    for(j in (i+1):length(raw.files)) {
      cat('\r', i, '-', j, '       ')  
      flush.console()
      
      ev.b <- ev %>%
        filter(PEP < pep.thresh) %>%
        filter(`Raw file`==raw.files[j]) %>%
        filter(`Mod. peptide ID` %in% ev.a$`Mod. peptide ID`) %>%
        select(c('Mod. peptide ID', 'Retention time')) %>%
        group_by(`Mod. peptide ID`) %>%
        summarise(`Retention time`=mean(`Retention time`))
      
      # require at least 5 common points to do this analysis
      if(nrow(ev.b) < 5) {
        cor.mat[i,j] <- NA
        next
      }
      
      b <- ev.b$`Retention time`
      a <- ev.a %>%
        filter(`Mod. peptide ID` %in% ev.b$`Mod. peptide ID`) %>%
        group_by(`Mod. peptide ID`) %>%
        summarise(`Retention time`=mean(`Retention time`)) %>%
        pull(`Retention time`)
      
      cor.mat[i,j] <- cor(a, b)
    }
  }
}

# only take upper triangle. don't want to destroy the matrix structure,
# so set the lower triangle to a dummy value so it can be removed later
cor.mat[lower.tri(cor.mat)] <- 100

# set diagnoal to 1
diag(cor.mat) <- 1

df <- melt(cor.mat)
df <- df %>%
  filter(!value == 100)

ggplot(df, aes(y=Var1, x=Var2, fill=value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", na.value='grey50',
                       midpoint = 0, limit = c(-1,1), space = "Lab",
                       name="Pearson\nCorrelation") +
  theme_minimal() +
  theme(axis.text.y = element_text(size=8),
        axis.text.x = element_text(angle = 45, vjust = 1,
                                   size = 8, hjust = 1)
        #axis.text.x = element_blank()
        ) +
  labs(x=NULL, y=NULL) +
  coord_fixed()
