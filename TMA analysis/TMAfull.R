library("tidyverse")
# Set the working directory for you TMA data and reference metadata
setwd("C:/Users/wp1g19/Documents/TMA analysis/R full analysis wd")
#Read in all data files into R environment
temp = list.files(pattern="*.csv")
for (i in 1:length(temp)) assign(temp[i], read.csv(temp[i]))

# Merge  counts/scores into a with reference sheet
CD3TMAreferenced <- merge(`CD3 TMA measurements.csv`, `Reference sheet.csv`, by=c("TMA.block","TMA.core"), )
CD8TMAreferenced <- merge(`CD8 TMA measurements.csv`, `Reference sheet.csv`, by=c("TMA.block","TMA.core"))
CD4TMAreferenced <- merge(`CD4 TMA measurements.csv`, `Reference sheet.csv`, by=c("TMA.block","TMA.core"))
FOXP3TMAreferenced <- merge(`FOXP3 TMA measurements.csv`, `Reference sheet.csv`, by=c("TMA.block","TMA.core"))
HLA_C2TMAreferenced <- merge(`HLA_Class_II TMA measurements.csv`, `Reference sheet.csv`, by=c("TMA.block","TMA.core"))
HLA_C2TMAreferenced_H_score <- merge(`HLA_Class_II H-Score TMA measurements.csv`, `Reference sheet.csv`, by=c("TMA.block","TMA.core"))
HLA_ABCTMAreferenced <- merge(`HLA_ABC TMA measurements.csv`, `Reference sheet.csv`, by=c("TMA.block","TMA.core"))
HLA_ABCTMAreferenced_H_score <- merge(`HLA_ABC H-Score TMA measurements.csv`, `Reference sheet.csv`, by=c("TMA.block","TMA.core"))
ERAP2TMAreferenced <- merge(`ERAP2 H-Score TMA measurements.csv`, `Reference sheet.csv`, by=c("TMA.block","TMA.core"))
HLA_ETMAreferenced <- merge(`HLA_E TMA measurements.csv`, `Reference sheet.csv`, by=c("TMA.block","TMA.core"))
HLA_ETMAreferenced_H_score <- merge(`HLA_E H-Score TMA measurements.csv`, `Reference sheet.csv`, by=c("TMA.block","TMA.core"))
TAP1TMAreferenced_H_score <- merge(`TAP1 H-Score TMA measurements.csv`, `Reference sheet.csv`, by=c("TMA.block","TMA.core"))

#Omit cases with low detections of overall cells/score, these will be the missing/bad cores
CD3TMAreferencedFilter <- CD3TMAreferenced %>% filter(Num.Detections > 50)
CD8TMAreferencedFilter <- CD8TMAreferenced %>% filter(Num.Detections > 50)
CD4TMAreferencedFilter <- CD4TMAreferenced %>% filter(Num.Detections > 50)
FOXP3TMAreferencedFilter <- FOXP3TMAreferenced %>% filter(Num.Detections > 50)
HLA_ABCTMAreferencedFilter <- HLA_ABCTMAreferenced #%>% filter(Total.ROI.area.µm.2 > 1000)
HLA_ETMAreferencedFilter <- HLA_ETMAreferenced #%>% filter(Total.ROI.area.µm.2 > 1000)
HLA_C2TMAreferencedFilter <- HLA_C2TMAreferenced #%>% filter(Total.ROI.area.µm.2 > 1000)
TAP1TMAreferenced_H_score <- TAP1TMAreferenced_H_score #%>% filter(Num.Detections > 100)

#Aggregate and calculate mean counts/score by histopathology number, relable columns
CD3countsbyhistNum <- group_by(CD3TMAreferencedFilter, Histopathology.Reference.Number) %>% summarize(m = mean(Num.Positive.per.mm.2))
colnames(CD3countsbyhistNum)  <- c("HistopathologyNumber", "CD3countpermm2")

CD8countsbyhistNum <- group_by(CD8TMAreferencedFilter, Histopathology.Reference.Number) %>% summarize(m = mean(Num.Positive.per.mm.2))
colnames(CD8countsbyhistNum)  <- c("HistopathologyNumber", "CD8countpermm2")

CD4countsbyhistNum <- group_by(CD4TMAreferencedFilter, Histopathology.Reference.Number) %>% summarize(m = mean(Num.Positive.per.mm.2))
colnames(CD4countsbyhistNum)  <- c("HistopathologyNumber", "CD4countpermm2")

FOXP3countsbyhistNum <- group_by(FOXP3TMAreferencedFilter, Histopathology.Reference.Number) %>% summarize(m = mean(Num.Positive.per.mm.2))
colnames(FOXP3countsbyhistNum)  <- c("HistopathologyNumber", "FOXP3countpermm2")

HLA_ABCscoresbyhistNum <- group_by(HLA_ABCTMAreferenced, Histopathology.Reference.Number) %>% summarize(m = mean(Positive...of.stained.pixels))
colnames(HLA_ABCscoresbyhistNum)  <- c("HistopathologyNumber", "HLA_ABCpercentagescore")

HLA_ABCHscoresbyhistNum <- group_by(HLA_ABCTMAreferenced_H_score, Histopathology.Reference.Number) %>% summarize(m = mean(H.score))
colnames(HLA_ABCHscoresbyhistNum)  <- c("HistopathologyNumber", "HLA_ABCHscore")

HLA_C2scoresbyhistNum <- group_by(HLA_C2TMAreferenced, Histopathology.Reference.Number) %>% summarize(m = mean(Positive...of.stained.pixels))
colnames(HLA_C2scoresbyhistNum)  <- c("HistopathologyNumber", "HLA_C2percentagescore")
HLA_C2HscoresbyhistNum <- group_by(HLA_C2TMAreferenced_H_score, Histopathology.Reference.Number) %>% summarize(m = mean(H.score))
colnames(HLA_C2HscoresbyhistNum)  <- c("HistopathologyNumber", "HLA_C2Hscore")


HLA_EscoresbyhistNum <- group_by(HLA_ETMAreferenced, Histopathology.Reference.Number) %>% summarize(m = mean(Positive...of.stained.pixels))
colnames(HLA_EscoresbyhistNum)  <- c("HistopathologyNumber", "HLA_Epercentagescore")
HLA_EHscoresbyhistNum <- group_by(HLA_ETMAreferenced_H_score, Histopathology.Reference.Number) %>% summarize(m = mean(H.score))
colnames(HLA_EHscoresbyhistNum)  <- c("HistopathologyNumber", "HLA_EHscore")

TAP1scoresbyhistNum<- group_by(TAP1TMAreferenced_H_score, Histopathology.Reference.Number) %>% summarize(m = mean(Positive...of.stained.pixels))
colnames(TAP1scoresbyhistNum)  <- c("HistopathologyNumber", "TAP1percentagescore")
TAP1TMAreferenced_H_score$H.score <- as.numeric(TAP1TMAreferenced_H_score$H.score)
TAP1HscoresbyhistNum <- group_by(TAP1TMAreferenced_H_score, Histopathology.Reference.Number) %>% summarize(m = mean(H.score))
colnames(TAP1HscoresbyhistNum)  <- c("HistopathologyNumber", "TAP1Hscore")

#Merge the data
df_list <- list(CD3countsbyhistNum, CD8countsbyhistNum ,CD4countsbyhistNum ,FOXP3countsbyhistNum, HLA_ABCscoresbyhistNum, HLA_ABCHscoresbyhistNum, HLA_C2scoresbyhistNum, HLA_C2HscoresbyhistNum, HLA_EscoresbyhistNum, HLA_EHscoresbyhistNum, TAP1scoresbyhistNum, TAP1HscoresbyhistNum)
Allmarkermerged <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list)
# remove the orient core
Allmarkermerged <- Allmarkermerged[!Allmarkermerged$HistopathologyNumber == "Orient",]
rownames(Allmarkermerged) <- Allmarkermerged[,1]
Allmarkermerged <- Allmarkermerged[,-1]
omitNAAllmarkerMerged <- na.omit(Allmarkermerged)
write.csv(omitNAAllmarkerMerged, file = "C:/Users/wp1g19/Documents/TMA analysis/Heatmapper/NAomitTMAanalysis for heatmapper.csv")

#Correlate the counts/scores
setwd("C:/Users/wp1g19/Documents/TMA analysis/R full analysis wd/images for thesis")
library("PerformanceAnalytics")
chart.Correlation(Allmarkermerged, histogram=TRUE, pch= "+", method = "pearson")
dev.copy(png,filename="Correlation.png");
dev.off ();
chart.Correlation(omitNAAllmarkerMerged, histogram=TRUE, pch= "+", method = "pearson")
dev.copy(png,filename="CorrelationOmitNA.png");
dev.off ();

#align Clinical data to the counts/scores
rownames(`TMA datasheet final referenced.csv`) <- `TMA datasheet final referenced.csv`[,2]
`TMA datasheet final referenced.csv` <- as.data.frame(`TMA datasheet final referenced.csv`)
Allmarkermerged$Histopathology.Reference.Number <- rownames(Allmarkermerged)

#save the merged counts/scores
write.csv(Allmarkermerged, file = "C:/Users/wp1g19/Documents/TMA analysis/R full analysis wd/Allmarkersmerged.csv")

#readback the csvs
Allmarkermerged <- read.csv(file = "C:/Users/wp1g19/Documents/TMA analysis/R full analysis wd/Allmarkersmerged.csv", row.names = 1)
Allmarkermerged %>% 
  mutate(across(where(is.character), str_remove_all, pattern = fixed(" ")))
TMACLINmerge <- merge(Allmarkermerged, `TMA datasheet final referenced.csv`, by = "Histopathology.Reference.Number")
#convert survival from days to months
TMACLINmerge <- TMACLINmerge %>% 
  mutate(OS.months = round(OS.survival.time/30.417, digit=0))

#Set WD to image file
setwd("C:/Users/wp1g19/Documents/TMA analysis/R full analysis wd/images for thesis")
#Overall survival with optimal cutoffs
library(survminer)
TMACLINmerge$HLA_EHscore
#make the optimal cut
HLA_EHscore.surv_count.cut <- surv_cutpoint(
  TMACLINmerge,
  time = "OS.months",
  event = "Survival.status",
  variables = "HLA_EHscore"
)

summary(HLA_EHscore.surv_count.cut)

#Plot the cutpoint
# palette = "npg" (nature publishing group), see ?ggpubr::ggpar
plot(HLA_EHscore.surv_count.cut, "HLA_EHscore", palette = "npg")
dev.copy(png,filename="HLA_EHscore.surv_count.cut.png");
dev.off ();

#categorise the cutpoint 
HLA_EHscore.surv_count.cat <- surv_categorize(HLA_EHscore.surv_count.cut) 

library(survival)
library(RTCGA)
fit <- survfit(Surv(OS.months, Survival.status) ~ HLA_EHscore,
               data = HLA_EHscore.surv_count.cat)
res <- pairwise_survdiff(Surv(OS.months, Survival.status) ~ HLA_EHscore,
                         data = HLA_EHscore.surv_count.cat, p.adjust.method = "fdr")
res
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,120), 
  xlab = ("Months"),# present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 20,    # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = TRUE # show bars instead of names in text annotations
  # in legend of risk table
)  
#save the image file
dev.copy(png,filename="HLA_EHscoreOSKm.png");
dev.off ();


#proportion of clinical features 
library(psych)
#Pull out the characteristics we want to tabulate
TMAdemoClin <- TMACLINmerge[, c("Age", "Sex", "pT", "pN", "pM", "LymphaticInvasion", "VascularInvasion", "Perineural_Invasion", "TumourGradingDifferentiation", "R0.American.", "MandardScore", "AliveorDead")]
# Catagorise missing data as "unknown"
TMAdemoClin2 <- as.data.frame(apply(TMAdemoClin, 2, function(y) (gsub("-", "Unknown", y))))

#define mean age and range
function(x) as.character(as.numeric(x))
TMAdemoClin2$Age <- as.numeric(as.character(TMAdemoClin2$Age))
library(psych)
Democlindescriptstats <- describe(TMAdemoClin2)

library(dplyr)
library(tidyr)
TMAdemoClin2
TMAdemoClin2 %>% 
  group_by(AliveorDead) %>% 
  tally()

write.csv(TMAdemoClin2, file = "C:/Users/wp1g19/Documents/TMA analysis/R full analysis wd/Clinicaldemographics.csv")




#Univariate CoxPH
covariates <- c("CD3countpermm2", "CD4countpermm2",  "CD8countpermm2", "FOXP3countpermm2")

CD3coxPH <- coxph(formula = Surv(OS.survival.time, Survival.status) ~ CD3countpermm2, data=TMACLINmerge)
CD8coxPH <- coxph(formula = Surv(OS.survival.time, Survival.status) ~ CD8countpermm2, data=TMACLINmerge)
CD4coxPH <- coxph(formula = Surv(OS.survival.time, Survival.status) ~ CD4countpermm2, data=TMACLINmerge)
FOXP3coxPH <- coxph(formula = Surv(OS.survival.time, Survival.status) ~ FOXP3countpermm2, data=TMACLINmerge)

covariate_names <- c(CD3density = "CD3countpermm2", 
                     CD8density = "CD8countpermm2",
                     CD4density = "CD4countpermm2",
                     FOXP3density = "FOXP3countpermm2")

map(vars(CD3countpermm2, CD8countpermm2, CD4countpermm2, FOXP3countpermm2), function(by)
{
  analyse_multivariate(TMACLINmerge,
                       vars(OS.survival.time, Survival.status),
                       covariates = list(by), # covariates expects a list
                       covariate_name_dict = covariate_names)
}) %>%
  forest_plot(factor_labeller = covariate_names,
              endpoint_labeller = c(time="OS"),
              orderer = ~order(HR),
              labels_displayed = c("endpoint", "factor", "n"),
              ggtheme = ggplot2::theme_bw(base_size = 10))

#Pull out the OS, censor and score for each stain

CD3extract <- TMACLINmerge[,c("Histopathology.Reference.Number","OS.months","Survival.status", "CD3countpermm2")]
colnames(CD3extract) <- c("Sample ID", "Survival time", "Survival event", "CD3countpermm2")

CD4extract <- TMACLINmerge[,c("Histopathology.Reference.Number","OS.months","Survival.status","CD4countpermm2")]
colnames(CD4extract) <- c("Sample ID", "Survival time", "Survival event", "CD4countpermm2")

CD8extract <- TMACLINmerge[,c("Histopathology.Reference.Number","OS.months","Survival.status","CD8countpermm2")]
colnames(CD8extract) <- c("Sample ID", "Survival time", "Survival event", "CD8countpermm2")

FOXP3extract <- TMACLINmerge[,c("Histopathology.Reference.Number","OS.months","Survival.status","FOXP3countpermm2")]
colnames(FOXP3extract) <- c("Sample ID", "Survival time", "Survival event", "FOXP3countpermm2")

HLA_ABCHscoreextract <- TMACLINmerge[,c("Histopathology.Reference.Number","OS.months","Survival.status", "HLA_ABCHscore")]
colnames(HLA_ABCHscoreextract) <- c("Sample ID", "Survival time", "Survival event", "HLA_ABCHscore")

HLA_ABCperextract <- TMACLINmerge[,c("Histopathology.Reference.Number","OS.months","Survival.status", "HLA_ABCpercentagescore")]
colnames(HLA_ABCperextract) <- c("Sample ID", "Survival time", "Survival event", "HLA_ABCpercentagescore")

HLA_EHscoreextract <- TMACLINmerge[,c("Histopathology.Reference.Number","OS.months","Survival.status","HLA_EHscore")]
colnames(HLA_EHscoreextract) <- c("Sample ID", "Survival time", "Survival event", "HLA_EHscore")

HLA_Eperextract <- TMACLINmerge[,c("Histopathology.Reference.Number","OS.months","Survival.status","HLA_Epercentagescore")]
colnames(HLA_Eperextract) <- c("Sample ID", "Survival time", "Survival event", "HLA_Epercentagescore")

HLA_C2Hscoreextract <- TMACLINmerge[,c("Histopathology.Reference.Number","OS.months","Survival.status","HLA_C2Hscore")]
colnames(HLA_C2Hscoreextract) <- c("Sample ID", "Survival time", "Survival event", "HLA_C2Hscore")

HLA_C2perextract <- TMACLINmerge[,c("Histopathology.Reference.Number","OS.months","Survival.status","HLA_C2percentagescore")]
colnames(HLA_C2perextract) <- c("Sample ID", "Survival time", "Survival event", "HLA_C2percentagescore")


TAP1Hscoreextract <- TMACLINmerge[,c("Histopathology.Reference.Number","OS.months","Survival.status","TAP1Hscore")]
colnames(TAP1Hscoreextract) <- c("Sample ID", "Survival time", "Survival event", "TAP1Hscore")

TAP1perextract <- TMACLINmerge[,c("Histopathology.Reference.Number","OS.months","Survival.status", "TAP1percentagescore")]
colnames(TAP1perextract) <- c("Sample ID", "Survival time", "Survival event", "TAP1percentagescore")

df_list <- list(CD3extract,
                CD4extract,
                CD8extract,
                FOXP3extract,
                HLA_ABCHscoreextract,
                HLA_ABCperextract,
                HLA_EHscoreextract,
                HLA_Eperextract,
                HLA_C2Hscoreextract,
                HLA_C2perextract,
                TAP1Hscoreextract,
                TAP1perextract)
Allmarkerextract <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list)



setwd("C:/Users/wp1g19/Documents/TMA analysis/R full analysis wd/extractedOS")
data_names <- c("CD3extract",
               "CD4extract",
               "CD8extract",
               "FOXP3extract",
               "HLA_ABCHscoreextract",
               "HLA_ABCperextract",
               "HLA_EHscoreextract",
               "HLA_Eperextract",
               "HLA_C2Hscoreextract",
               "HLA_C2perextract",
               "TAP1Hscoreextract",
               "TAP1perextract",
               "Allmarkerextract")


for(i in 1:length(data_names)) {                              # Head of for-loop
  write.csv(get(data_names[i]),                              # Write CSV files to folder
             paste0("C:/Users/wp1g19/Documents/TMA analysis/R full analysis wd/extractedOS/",
                    data_names[i],
                    ".csv"),
             row.names = FALSE)
}

