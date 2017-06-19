library(scales)
library(readr)

evidence <- read_tsv("~/Google Drive/Ali_RT_bayesian/dat/evidence.txt")

hist(evidence[["Retention time"]], breaks=50, probability=TRUE)
abline(v=mean(evidence[evidence[["Peptide ID"]] == 208, ][["Retention time"]]),
       lwd=3, col="blue")
abline(v=mean(evidence[evidence[["Peptide ID"]] == 290, ][["Retention time"]]),
       lwd=3, col="red")
abline(v=mean(evidence[evidence[["Peptide ID"]] == 196, ][["Retention time"]]),
       lwd=3, col="green")


which.max(table(evidence[["Leading razor protein"]]))

sub <- evidence[evidence[["Leading razor protein"]] == "sp|P36578|RL4_HUMAN", ]

sub[["Peptide ID"]]

######################
## Peptide Exploration
######################

peptideCounts <- sort(table(evidence["Peptide ID"]), decreasing=TRUE)
hist(log10(as.numeric(peptideCounts)))
head(peptideCounts, n=100)

## to look at
## 592348, 572738

## explore some peptides
pep <- 592348
fdr <- pmin(evidence$PEP[evidence["Peptide ID"] == pep], 1)
print(fdr)

## hist(evidence[["Retention time"]], breaks=50, probability=TRUE, ylim=c(0, 0.01))
rts <- evidence[evidence[["Peptide ID"]] == pep, ][["Retention time"]]
hist(evidence[["Retention time"]], breaks=50, probability=TRUE)
abline(rts,col=alpha("blue", alpha=(1-fdr)))

######################
## Protein Exploration
######################

proteinCounts <- sort(table(evidence["Leading razor protein"]), decreasing=TRUE)
## for(protein in names(proteinCounts)) {
protein <- "sp|P19338|NUCL_HUMAN"
proteinEvidence <- evidence[evidence["Leading razor protein"] == protein, ]
table(proteinEvidence[["Peptide ID"]])

## 474, 480 , 467437
i1 <- proteinEvidence[proteinEvidence["Peptide ID"] == 474, c("Raw file", "Intensity", "PEP")]
i2 <- proteinEvidence[proteinEvidence["Peptide ID"] == 418167, c("Raw file", "Intensity", "PEP")]

indices <- intersect(i1[["Raw file"]], i2[["Raw file"]])
res1 <- log(i1[match(indices, i1[["Raw file"]]), "Intensity"])
res2 <- log(i2[match(indices, i2[["Raw file"]]), "Intensity"])

pep1 <- pmin(i1[match(indices, i1[["Raw file"]]), ]$PEP, 1)
pep2 <- pmin(i2[match(indices, i2[["Raw file"]]), ]$PEP, 1)
sort(pep1)
sort(pep2)
cor.test(pep1, pep2)

cor.test(res1$Intensity[-1], res2$Intensity[-1], pch=19)
plot(res1$Intensity, res2$Intensity, pch=19)
plot(res1$Intensity, res2$Intensity, pch=19, col=alpha("black", alpha=abs(pep1-pep2)))


##}

## for( rf in unique(evidence[["Raw file"]])) {
##     print(summary(log(evidence$Intensity[evidence[["Raw file"]] == rf]), na.rm=TRUE))
## }
