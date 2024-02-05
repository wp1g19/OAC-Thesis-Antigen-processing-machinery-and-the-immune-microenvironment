#TMA EDA
library(tidyverse)
library(psych)
#read in csv file
cd3data <- read.csv(file = "C:/Users/wp1g19/Documents/TMA analysis/CD3 TMA measurements.csv", header = TRUE)

cd4data <- read.csv(file = "C:/Users/wp1g19/Documents/TMA analysis/CD4 TMA measurements.csv", header = TRUE)

cd8data <- read.csv(file = "C:/Users/wp1g19/Documents/TMA analysis/CD8 TMA measurements.csv", header = TRUE)

foxp3data <- read.csv(file = "C:/Users/wp1g19/Documents/TMA analysis/FOXP3 TMA measurements.csv", header = TRUE)

HLA_ABC <- read.csv(file = "C:/Users/wp1g19/Documents/TMA analysis/HLA_ABC TMA measurements.csv", header = TRUE)

HLA_E <- read.csv(file = "C:/Users/wp1g19/Documents/TMA analysis/HLA_E TMA measurements.csv", header = TRUE)

TAP1 <- read.csv(file = "C:/Users/wp1g19/Documents/TMA analysis/TAP1 TMA measurements.csv", header = TRUE)

HLA_Class_II <- read.csv(file = "C:/Users/wp1g19/Documents/TMA analysis/HLA_Class_II TMA measurements.csv", header = TRUE)
#filtering and quality control
# remove samples with zero detections
cd3data = subset(cd3data, cd3data$Num.Positive.per.mm.2 > 0)
cd3data = subset(cd3data, cd3data$Num.Detections > 500)

cd4data = subset(cd4data, cd4data$Num.Positive.per.mm.2 > 0)
cd4data = subset(cd4data, cd4data$Num.Detections > 500)

cd8data = subset(cd8data, cd8data$Num.Positive.per.mm.2 > 0)
cd8data = subset(cd8data, cd8data$Num.Detections > 500)

foxp3data = subset(foxp3data, foxp3data$Num.Positive.per.mm.2 > 0)
foxp3data = subset(foxp3data, foxp3data$Num.Detections > 500)

HLA_ABC = subset(HLA_ABC, HLA_ABC$Positive.pixel.area.Âµm.2 > 1000)

HLA_E = subset(HLA_E, HLA_E$Positive.pixel.area.Âµm.2 > 1000)

HLA_Class_II = subset(HLA_Class_II, HLA_Class_II$Positive.pixel.area.Âµm.2 > 1000)

TAP1 = subset(TAP1, TAP1$Positive.pixel.area.µm.2..d.4..s.2..tN.0.1..tP.0.3. > 100)

#description of data

cd3datasummary <- as.data.frame(summary(cd3data))

cd4datasummary <- as.data.frame(summary(cd4data))

cd8datasummary <- as.data.frame(summary(cd8data))

foxp3datasummary <- as.data.frame(summary(foxp3data))

HLA_ABCdatasummary <- as.data.frame(summary(HLA_ABC))

HLA_Edatasummary <- as.data.frame(summary(HLA_E))

HLA_C2datasummary <- as.data.frame(summary(HLA_Class_II))

cd3datadescript <- as.data.frame(describe(cd3data))

cd4datadescript <- as.data.frame(describe(cd4data))

cd8datadescript <- as.data.frame(describe(cd8data))

foxp3datadescript <- as.data.frame(describe(foxp3data))

HLA_ABCdatadescript <- as.data.frame(describe(HLA_ABC))

HLA_Edatadescript <- as.data.frame(describe(HLA_E))

HLA_C2datasummary <- as.data.frame(summary(HLA_Class_II))

TAP1datasummary <- as.data.frame(summary(TAP1))

#Histogram the distribution of immune cell counts

CD3countpermm2 <- cd3data$Num.Positive.per.mm.2
hist(CD3countpermm2)

CD4countpermm2 <- cd4data$Num.Positive.per.mm.2
hist(CD4countpermm2)

CD8countpermm2 <- cd8data$Num.Positive.per.mm.2
hist(CD8countpermm2)

foxp3countpermm2 <- foxp3data$Num.Positive.per.mm.2
hist(foxp3countpermm2)

HLA_ABC_perc_score <- HLA_ABC$Positive...of.stained.pixels
hist(HLA_ABC_perc_score)

HLA_E_perc_score <- HLA_E$Positive...of.stained.pixels
hist(HLA_E_perc_score)

HLA_C2_perc_score <- HLA_Class_II$Positive...of.stained.pixels
hist(HLA_C2_perc_score)

TAP1_perc_score <- TAP1$Positive...of.stained.pixels..d.4..s.2..tN.0.1..tP.0.3.
hist(TAP1_perc_score)
#summary and description of counts/scores

CD3countpermm2df <- as.data.frame(CD3countpermm2)

CD4countpermm2df <- as.data.frame(CD4countpermm2)

CD8countpermm2df <- as.data.frame(CD8countpermm2)

foxp3countpermm2df <- as.data.frame(foxp3countpermm2)

HLA_ABC_perc_scoredf <- as.data.frame(HLA_ABC_perc_score)

HLA_E_perc_scoredf <- as.data.frame(HLA_E_perc_score)

HLA_C2_perc_scoredf <- as.data.frame(HLA_C2_perc_score)

TAP1_perc_scoredf <- as.data.frame(TAP1_perc_score)

cd3dcountatasummary <- as.data.frame(summary(CD3countpermm2df))

cd4datacountsummary <- as.data.frame(summary(CD4countpermm2df))

cd8countdatasummary <- as.data.frame(summary(CD8countpermm2df))

foxp3countdatasummary <- as.data.frame(summary(foxp3countpermm2df))

HLA_ABCpercscoredatasummary <- as.data.frame(summary(HLA_ABC_perc_scoredf))

HLA_Eperccoredatasummary <- as.data.frame(summary(HLA_E_perc_scoredf))

HLA_C2perccoredatasummary <- as.data.frame(summary(HLA_C2_perc_scoredf))

cd3countdatadescript <- as.data.frame(describe(CD3countpermm2df))

cd4countdatadescript <- as.data.frame(describe(CD4countpermm2df))

cd8countdatadescript <- as.data.frame(describe(CD8countpermm2df))

foxp3countdatadescript <- as.data.frame(describe(foxp3countpermm2df))

HLA_ABCpercscoredatadescript <- as.data.frame(describe(HLA_ABC_perc_scoredf))

HLA_Epercscoredatadescript <- as.data.frame(describe(HLA_E_perc_scoredf))

HLA_C2percscoredatadescript <- as.data.frame(describe(HLA_C2_perc_scoredf))

TAP1percscoredatadescript <- as.data.frame(describe(TAP1_perc_scoredf))

#histogram with standard curve

hist(CD3countpermm2, freq = FALSE, main = "Density curve", ylim = c(0,0.01), xlim = c(0,2500), breaks = 381)
lines(density(CD3countpermm2), lwd = 2, col = 'red')

hist(CD4countpermm2, freq = FALSE, main = "Density curve", ylim = c(0,0.01), xlim = c(0,2500), breaks = 402)
lines(density(CD4countpermm2), lwd = 2, col = 'red')

hist(CD8countpermm2, freq = FALSE, main = "Density curve", ylim = c(0,0.01), xlim = c(0,2500), breaks = 431)
lines(density(CD8countpermm2), lwd = 2, col = 'red')

hist(foxp3countpermm2, freq = FALSE, main = "Density curve", breaks = 392)
lines(density(foxp3countpermm2), lwd = 2, col = 'red')

hist(HLA_ABC_perc_score, freq = FALSE, main = "Density curve",ylim = c(0,0.10), breaks = 402)
lines(density(HLA_ABC_perc_score), lwd = 2, col = 'red')

hist(HLA_E_perc_score, freq = FALSE, main = "Density curve", breaks = 301, ylim = c(0,0.30), xlim = c(0,100))
lines(density(HLA_E_perc_score), lwd = 2, col = 'red')

hist(HLA_C2_perc_score, freq = FALSE, main = "Density curve", breaks = 384, xlim = c(0,100), ylim = c(0,0.2), probability = TRUE)
lines(density(HLA_C2_perc_score), lwd = 2, col = 'red')

hist(TAP1_perc_score, freq = FALSE, main = "Density curve", breaks = 384, xlim = c(0,100), ylim = c(0,0.2), probability = TRUE)
lines(density(TAP1_perc_score), lwd = 2, col = 'red')

#Histogram of H scores
HLA_ABC <- read.csv(file = "C:/Users/wp1g19/Documents/TMA analysis/HLA_ABC H-Score TMA measurements.csv", header = TRUE)

HLA_E <- read.csv(file = "C:/Users/wp1g19/Documents/TMA analysis/HLA_E H-Score TMA measurements.csv", header = TRUE)

HLA_Class_II <- read.csv(file = "C:/Users/wp1g19/Documents/TMA analysis/HLA_Class_II H-Score TMA measurements.csv", header = TRUE)

TAP1 <- TAP1

HLA_ABC = subset(HLA_ABC, HLA_ABC$Positive.pixel.area.Âµm.2 > 1000)

HLA_E = subset(HLA_E, HLA_E$Positive.pixel.area.Âµm.2 > 1000)

HLA_Class_II = subset(HLA_Class_II, HLA_Class_II$Positive.pixel.area.Âµm.2 > 1000)

TAP1 = subset(TAP1, HLA_Class_II$H.Score�µm.2 > 1)

TAP1$H.score <- as.numeric(TAP1$H.score)

HLA_ABC_H_score <- HLA_ABC$H.score

HLA_E_H_score <- HLA_E$H.score

HLA_C2_H_score <- HLA_Class_II$H.Score

TAP1_H_score <- TAP1$H.score

hist(HLA_ABC_H_score, freq = FALSE, main = "Density curve",ylim = c(0,0.10), breaks = 402)
lines(density(HLA_ABC_H_score), lwd = 2, col = 'red')

hist(HLA_E_H_score, freq = FALSE, main = "Density curve", breaks = 301, ylim = c(0,0.30), xlim = c(0,100))
lines(density(HLA_E_H_score), lwd = 2, col = 'red')

hist(HLA_C2_H_score, freq = FALSE, main = "Density curve", breaks = 384, xlim = c(0,100), ylim = c(0,0.2), probability = TRUE)
lines(density(HLA_C2_H_score), lwd = 2, col = 'red')

hist(TAP1_H_score, freq = FALSE, main = "Density curve", breaks = 384, xlim = c(0,100), ylim = c(0,0.2), probability = TRUE)
lines(density(TAP1_H_score), lwd = 2, col = 'red')

#summary of H-scores

HLA_ABC_h_scoredf <- as.data.frame(HLA_ABC$H.score)

HLA_E_h_scoredf <- as.data.frame(HLA_E$H.score)

HLA_C2_h_scoredf <- as.data.frame(HLA_Class_II$H.Score)

TAP1_h_scoredf <- as.data.frame(TAP1_H_score)

HLA_ABCHscoredatasummary <- as.data.frame(summary(HLA_ABC_h_scoredf))

HLA_EHcoredatasummary <- as.data.frame(summary(HLA_E_h_scoredf))

HLA_C2Hcoredatasummary <- as.data.frame(summary(HLA_C2_h_scoredf))

HLA_ABCHscoredatadescript <- as.data.frame(describe(HLA_ABC_h_scoredf))

HLA_EHscoredatadescript <- as.data.frame(describe(HLA_E_h_scoredf))

TAP1Hscoredatadescript <- as.data.frame(describe(TAP1_h_scoredf))

HLA_C2Hscoredatadescript <- as.data.frame(describe(HLA_C2_h_scoredf))

summary(TAP1_perc_scoredf)
summary(TAP1_h_scoredf)
