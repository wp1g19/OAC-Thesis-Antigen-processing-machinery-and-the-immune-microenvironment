#high HLA-ABC only
library(tidyverse)
abchighonly <- read.csv(file = "C:/Users/wp1g19/Documents/TMA analysis/R full analysis wd/Maximal concordant survival analysis/OS_APM_Immunepartners/CutpointsHLA_ABCHighOnly.csv")
PMDATA <- COXphoptimalcutoffsClinMerge[c("pM", "Histopathology.Reference.Number")]
library(finalfit)

abchighonly <- merge(abchighonly, PMDATA, by = "Histopathology.Reference.Number")
dependent_os  <- "Surv(OS.months, Survival.status)"
abchighonly$HLA_ABC_CD8
explanatory   <- c(    "Age",
                       "Sex",
                       "pT",
                       "pN",
                       "pM",
                       "HLA_ABC_CD8")


abchighonly %>% 
  finalfit(dependent_os, explanatory, metrics = TRUE)

allcuts <- read.csv(file = "C:/Users/wp1g19/Documents/TMA analysis/R full analysis wd/Maximal concordant survival analysis/OS_APM_Immunepartners/OptimalcutpointsALL.csv")
allcutsHLA_Ehigh <- allcuts[allcuts$HLA_EHscore == "high",]
allcutsHLA_Ehigh <- allcutsHLA_Ehigh %>% drop_na(HLA_EHscore)

allcutsHLA_Ehigh %>% 
  finalfit(dependent_os, explanatory, metrics = TRUE)

allcuts <- read.csv(file = "C:/Users/wp1g19/Documents/TMA analysis/R full analysis wd/Maximal concordant survival analysis/OS_APM_Immunepartners/OptimalcutpointsALL.csv")
allcutsTAP1high <- allcuts[allcuts$TAP1Hscore == "high",]
allcutsTAP1high <- allcutsTAP1high %>% drop_na(TAP1Hscore)

allcutsTAP1high %>% 
  finalfit(dependent_os, explanatory, metrics = TRUE)
