library("tidyverse")
# Set the working directory for you TMA data and reference metadata
setwd("C:/Users/wp1g19/Documents/TMA analysis/R full analysis wd/Measurements")
#Read in all data files into R environment
temp = list.files(pattern="*.csv")
for (i in 1:length(temp)) assign(temp[i], read.csv(temp[i]))

# Merge  counts/scores into a with reference sheet
CD3TMAreferenced <- merge(`CD3 TMA measurements.csv`, `Reference sheet.csv`, by=c("TMA.block","TMA.core"))
CD8TMAreferenced <- merge(`CD8 TMA measurements.csv`, `Reference sheet.csv`, by=c("TMA.block","TMA.core"))
CD4TMAreferenced <- merge(`CD4 TMA measurements.csv`, `Reference sheet.csv`, by=c("TMA.block","TMA.core"))
FOXP3TMAreferenced <- merge(`FOXP3 TMA measurements.csv`, `Reference sheet.csv`, by=c("TMA.block","TMA.core"))
HLA_C2TMAreferenced <- merge(`HLA_Class_II TMA measurements.csv`, `Reference sheet.csv`, by=c("TMA.block","TMA.core"))
HLA_C2TMAreferenced_H_score <- merge(`HLA_Class_II H-Score TMA measurements.csv`, `Reference sheet.csv`, by=c("TMA.block","TMA.core"))
HLA_ABCTMAreferenced <- merge(`HLA_ABC TMA measurements.csv`, `Reference sheet.csv`, by=c("TMA.block","TMA.core"))
HLA_ABCTMAreferenced_H_score <- merge(`HLA_ABC H-Score TMA measurements.csv`, `Reference sheet.csv`, by=c("TMA.block","TMA.core"))
HLA_ETMAreferenced <- merge(`HLA_E TMA measurements.csv`, `Reference sheet.csv`, by=c("TMA.block","TMA.core"))
HLA_ETMAreferenced_H_score <- merge(`HLA_E H-Score TMA measurements.csv`, `Reference sheet.csv`, by=c("TMA.block","TMA.core"))
TAP1TMAreferenced <- merge(`TAP1 TMA measurements.csv`, `Reference sheet.csv`, by=c("TMA.block","TMA.core"))
TAP1TMAreferenced_H_score <- merge(`TAP1 H-Score TMA measurements.csv`, `Reference sheet.csv`, by=c("TMA.block","TMA.core"))
CSDE1TMAreferenced_H_score <- merge(`CSDE1 H score measurements.csv`, `Reference sheet.csv`, by=c("TMA.block","TMA.core"))
CSDE1TMAreferenced <- merge(`CSDE1 measurements.csv`, `Reference sheet.csv`, by=c("TMA.block","TMA.core"))
#normalise csde1 to tumour content
classifiercsde1 <- read.csv(file = "C:/Users/wp1g19/Documents/TMA analysis/R full analysis wd/Measurements/CSDE1 superpixel measurements.csv", header = TRUE)

HLACSDE1Tumour <- merge(classifiercsde1, `Reference sheet.csv`, by=c("TMA.block","TMA.core"))
HLACSDE1Tumour$Percentage.Tumour <- as.numeric(HLACSDE1Tumour$Percentage.Tumour)
HLACSDE1Tumour <- group_by(HLACSDE1Tumour, Histopathology.Reference.Number) %>% summarize(m = mean(Percentage.Tumour))
colnames(HLACSDE1Tumour)  <- c("Histopathology.Reference.Number", "Percentage.Tumour")
HLACSDE1Tumour <- merge(HLACSDE1Tumour, CSDE1TMAreferenced_H_score, by = "Histopathology.Reference.Number")
HLACSDE1Tumour$H.score
  library(ggpubr)
ggplot(HLACSDE1Tumour, aes(x=Percentage.Tumour, y=H.score)) + 
  geom_point() +stat_smooth(method = "lm", col = "red") +
  stat_cor(label.x = 3, label.y = 250) 
write.csv(HLACSDE1Tumour, file = "C:/Users/wp1g19/Documents/TMA analysis/R full analysis wd/CSDE1HscoreTumourpercentage.csv")

#Omit cases with low detections of overall cells/score, these will be the missing/bad cores
CD3TMAreferencedFilter <- CD3TMAreferenced #%>% filter(Num.Detections > 50)
CD8TMAreferencedFilter <- CD8TMAreferenced #%>% filter(Num.Detections > 50)
CD4TMAreferencedFilter <- CD4TMAreferenced #%>% filter(Num.Detections > 50)
FOXP3TMAreferencedFilter <- FOXP3TMAreferenced #%>% filter(Num.Detections > 50)
HLA_ABCTMAreferencedFilter <- HLA_ABCTMAreferenced #%>% filter(Total.ROI.area.µm.2 > 1000)
HLA_ETMAreferencedFilter <- HLA_ETMAreferenced #%>% filter(Total.ROI.area.µm.2 > 1000)
HLA_ETMAreferencedFilter_H_score <- HLA_ETMAreferenced_H_score #%>% filter(Total.ROI.area.µm.2 > 1000)
HLA_C2TMAreferencedFilter <- HLA_C2TMAreferenced #%>% filter(Total.ROI.area.µm.2 > 1000)
HLA_C2TMAreferencedFilter_H_score <- HLA_C2TMAreferenced_H_score #%>% filter(Total.ROI.area.µm.2 > 1000)
TAP1TMAreferencedFilter_H_score <- TAP1TMAreferenced_H_score #%>% filter(Num.Detections > 1000)
TAP1TMAreferencedFilter <- TAP1TMAreferenced #%>% filter(Num.Detections > 1000)
CSDE1TMAreferencedFilter_H_score <- CSDE1TMAreferenced_H_score #%>% filter(Num.Detections > 1000)
CSDE1TMAreferencedFilter <- CSDE1TMAreferenced[CSDE1TMAreferenced$Histopathology.Reference.Number %in% CSDE1TMAreferencedFilter_H_score$Histopathology.Reference.Number,]
#Addition filters
#Quantile exclusion
Q <- quantile(CD3TMAreferencedFilter$Num.Positive.per.mm.2, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(CD3TMAreferencedFilter$Num.Positive.per.mm.2)
up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range???
CD3TMAreferencedFilter<- subset(CD3TMAreferencedFilter, CD3TMAreferencedFilter$Num.Positive.per.mm.2 > (Q[1] - 1.5*iqr) & CD3TMAreferencedFilter$Num.Positive.per.mm.2 < (Q[2]+1.5*iqr))

Q <- quantile(CD4TMAreferencedFilter$Num.Positive.per.mm.2, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(CD4TMAreferencedFilter$Num.Positive.per.mm.2)
up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range
CD4TMAreferencedFilter<- subset(CD4TMAreferencedFilter, CD4TMAreferencedFilter$Num.Positive.per.mm.2 > (Q[1] - 1.5*iqr) & CD4TMAreferencedFilter$Num.Positive.per.mm.2 < (Q[2]+1.5*iqr))


Q <- quantile(CD8TMAreferencedFilter$Num.Positive.per.mm.2, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(CD8TMAreferencedFilter$Num.Positive.per.mm.2)
up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range
CD8TMAreferencedFilter<- subset(CD8TMAreferencedFilter, CD8TMAreferencedFilter$Num.Positive.per.mm.2 > (Q[1] - 1.5*iqr) & CD8TMAreferencedFilter$Num.Positive.per.mm.2 < (Q[2]+1.5*iqr))

Q <- quantile(FOXP3TMAreferencedFilter$Num.Positive.per.mm.2, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(FOXP3TMAreferencedFilter$Num.Positive.per.mm.2)
up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range
FOXP3TMAreferencedFilter<- subset(FOXP3TMAreferencedFilter, FOXP3TMAreferencedFilter$Num.Positive.per.mm.2 > (Q[1] - 1.5*iqr) & FOXP3TMAreferencedFilter$Num.Positive.per.mm.2 < (Q[2]+1.5*iqr))

Q <- quantile(HLA_ABCTMAreferencedFilter$Positive...of.stained.pixels, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(HLA_ABCTMAreferencedFilter$Positive...of.stained.pixels)
up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range
HLA_ABCTMAreferencedFilter<- subset(HLA_ABCTMAreferencedFilter, HLA_ABCTMAreferencedFilter$Positive...of.stained.pixels > (Q[1] - 1.5*iqr) & HLA_ABCTMAreferencedFilter$Positive...of.stained.pixels < (Q[2]+1.5*iqr))

Q <- quantile(HLA_ETMAreferencedFilter$Positive...of.stained.pixels, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(HLA_ETMAreferencedFilter$Positive...of.stained.pixels)
up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range
HLA_ETMAreferencedFilter<- subset(HLA_ETMAreferencedFilter, HLA_ETMAreferencedFilter$Positive...of.stained.pixels > (Q[1] - 1.5*iqr) & HLA_ETMAreferencedFilter$Positive...of.stained.pixels < (Q[2]+1.5*iqr))

Q <- quantile(HLA_C2TMAreferencedFilter$Positive...of.stained.pixels, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(HLA_C2TMAreferencedFilter$Positive...of.stained.pixels)
up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range
HLA_C2TMAreferencedFilter<- subset(HLA_C2TMAreferencedFilter, HLA_C2TMAreferencedFilter$Positive...of.stained.pixels > (Q[1] - 1.5*iqr) & HLA_C2TMAreferencedFilter$Positive...of.stained.pixels < (Q[2]+1.5*iqr))

Q <- quantile(HLA_C2TMAreferencedFilter_H_score$H.score, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(HLA_C2TMAreferencedFilter_H_score$H.score)
up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range
HLA_C2TMAreferencedFilter_H_score<- subset(HLA_C2TMAreferencedFilter_H_score, HLA_C2TMAreferencedFilter_H_score$H.score > (Q[1] - 1.5*iqr) & HLA_C2TMAreferencedFilter_H_score$H.score < (Q[2]+1.5*iqr))

Q <- quantile(TAP1TMAreferenced_H_score$H.score, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(TAP1TMAreferenced_H_score$H.score)
up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range
TAP1TMAreferenced_H_score<- subset(TAP1TMAreferenced_H_score, TAP1TMAreferenced_H_score$H.score > (Q[1] - 1.5*iqr) & TAP1TMAreferenced_H_score$H.score < (Q[2]+1.5*iqr))

Q <- quantile(TAP1TMAreferencedFilter_H_score$Positive...of.stained.pixels, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(TAP1TMAreferencedFilter_H_score$Positive...of.stained.pixels)
up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range
TAP1TMAreferencedFilter_H_score<- subset(TAP1TMAreferencedFilter_H_score, TAP1TMAreferencedFilter_H_score$Positive...of.stained.pixels > (Q[1] - 1.5*iqr) & TAP1TMAreferencedFilter_H_score$Positive...of.stained.pixels < (Q[2]+1.5*iqr))

Q <- quantile(CSDE1TMAreferenced_H_score$H.score, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(CSDE1TMAreferenced_H_score$H.score)
up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range
CSDE1TMAreferencedFilter_H_score<- subset(CSDE1TMAreferencedFilter_H_score, CSDE1TMAreferencedFilter_H_score$H.score > (Q[1] - 1.5*iqr) & CSDE1TMAreferencedFilter_H_score$H.score < (Q[2]+1.5*iqr))

Q <- quantile(CSDE1TMAreferencedFilter$Positive...of.stained.pixels..d.4..s.2..tN.0.1..tP.0.3., probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(CSDE1TMAreferencedFilter$Positive...of.stained.pixels..d.4..s.2..tN.0.1..tP.0.3.)
up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range
CSDE1TMAreferencedFilter<- subset(CSDE1TMAreferencedFilter, CSDE1TMAreferencedFilter$Positive...of.stained.pixels..d.4..s.2..tN.0.1..tP.0.3. > (Q[1] - 1.5*iqr) & CSDE1TMAreferencedFilter$Positive...of.stained.pixels..d.4..s.2..tN.0.1..tP.0.3. < (Q[2]+1.5*iqr))

######################################################

#Aggregate and calculate mean counts/score by histopathology number, relable columns
CD3countsbyhistNum <- group_by(CD3TMAreferencedFilter, Histopathology.Reference.Number) %>% summarize(m = mean(Num.Positive.per.mm.2))
colnames(CD3countsbyhistNum)  <- c("HistopathologyNumber", "CD3countpermm2")

CD8countsbyhistNum <- group_by(CD8TMAreferencedFilter, Histopathology.Reference.Number) %>% summarize(m = mean(Num.Positive.per.mm.2))
colnames(CD8countsbyhistNum)  <- c("HistopathologyNumber", "CD8countpermm2")

CD4countsbyhistNum <- group_by(CD4TMAreferencedFilter, Histopathology.Reference.Number) %>% summarize(m = mean(Num.Positive.per.mm.2))
colnames(CD4countsbyhistNum)  <- c("HistopathologyNumber", "CD4countpermm2")

FOXP3countsbyhistNum <- group_by(FOXP3TMAreferencedFilter, Histopathology.Reference.Number) %>% summarize(m = mean(Num.Positive.per.mm.2))
colnames(FOXP3countsbyhistNum)  <- c("HistopathologyNumber", "FOXP3countpermm2")

HLA_ABCPercentagescoresbyhistNum <- group_by(HLA_ABCTMAreferenced, Histopathology.Reference.Number) %>% summarize(m = mean(Positive...of.stained.pixels))
colnames(HLA_ABCPercentagescoresbyhistNum)  <- c("HistopathologyNumber", "HLA_ABCpercentagescore")

HLA_ABCHscoresbyhistNum <- group_by(HLA_ABCTMAreferenced_H_score, Histopathology.Reference.Number) %>% summarize(m = mean(H.score))
colnames(HLA_ABCHscoresbyhistNum)  <- c("HistopathologyNumber", "HLA_ABCHscore")

HLA_C2TMAreferenced$Positive...of.stained.pixels <- as.numeric(HLA_C2TMAreferenced$Positive...of.stained.pixels)
HLA_C2PercentagescoresbyhistNum <- group_by(HLA_C2TMAreferenced, Histopathology.Reference.Number) %>% summarize(m = mean(Positive...of.stained.pixels))
colnames(HLA_C2PercentagescoresbyhistNum)  <- c("HistopathologyNumber", "HLA_C2 positive %")

HLA_C2scoresbyhistNum <- group_by(HLA_C2TMAreferenced, Histopathology.Reference.Number) %>% summarize(m = mean(Positive...of.stained.pixels))
colnames(HLA_C2scoresbyhistNum)  <- c("HistopathologyNumber", "HLA_C2percentagescore")
HLA_C2HscoresbyhistNum <- group_by(HLA_C2TMAreferenced_H_score, Histopathology.Reference.Number) %>% summarize(m = mean(H.score))
colnames(HLA_C2HscoresbyhistNum)  <- c("HistopathologyNumber", "HLA_C2Hscore")

HLA_ETMAreferencedFilter$Positive...of.stained.pixels <- as.numeric(HLA_ETMAreferencedFilter$Positive...of.stained.pixels)
HLA_EPercentagescoresbyhistNum <- group_by(HLA_ETMAreferencedFilter, Histopathology.Reference.Number) %>% summarize(m = mean(Positive...of.stained.pixels))
colnames(HLA_EPercentagescoresbyhistNum)  <- c("HistopathologyNumber", "HLA_E positive %")

HLA_EscoresbyhistNum <- group_by(HLA_ETMAreferenced, Histopathology.Reference.Number) %>% summarize(m = mean(Positive...of.stained.pixels))
colnames(HLA_EscoresbyhistNum)  <- c("HistopathologyNumber", "HLA_Epercentagescore")
HLA_EHscoresbyhistNum <- group_by(HLA_ETMAreferenced_H_score, Histopathology.Reference.Number) %>% summarize(m = mean(H.score))
colnames(HLA_EHscoresbyhistNum)  <- c("HistopathologyNumber", "HLA_EHscore")


TAP1TMAreferenced_H_score$H.score <- as.numeric(TAP1TMAreferenced_H_score$H.score)
TAP1HscoresbyhistNum <- group_by(TAP1TMAreferenced_H_score, Histopathology.Reference.Number) %>% summarize(m = mean(H.score))
colnames(TAP1HscoresbyhistNum)  <- c("HistopathologyNumber", "TAP1Hscore")

TAP1TMAreferenced$Positive...of.stained.pixels <- as.numeric(TAP1TMAreferenced$Positive...of.stained.pixels)
TAP1PercentagescoresbyhistNum <- group_by(TAP1TMAreferenced, Histopathology.Reference.Number) %>% summarize(m = mean(Positive...of.stained.pixels))
colnames(TAP1PercentagescoresbyhistNum)  <- c("HistopathologyNumber", "TAP1 positive %")


CSDE1TMAreferencedFilter_H_score$H.score <- as.numeric(CSDE1TMAreferencedFilter_H_score$H.score)
CSDE1scoresbyhistNum <- group_by(CSDE1TMAreferencedFilter_H_score, Histopathology.Reference.Number) %>% summarize(m = mean(H.score))
colnames(CSDE1scoresbyhistNum)  <- c("HistopathologyNumber", "CSDE1Hscore")

CSDE1TMAreferencedFilter$Positive...of.stained.pixels..d.4..s.2..tN.0.1..tP.0.3. <- as.numeric(CSDE1TMAreferencedFilter$Positive...of.stained.pixels..d.4..s.2..tN.0.1..tP.0.3.)
CSDE1PercentagescoresbyhistNum <- group_by(CSDE1TMAreferencedFilter, Histopathology.Reference.Number) %>% summarize(m = mean(Positive...of.stained.pixels..d.4..s.2..tN.0.1..tP.0.3.))
colnames(CSDE1PercentagescoresbyhistNum)  <- c("HistopathologyNumber", "CSDE1 positive %")
#Merge the data
df_list <- list(CD3countsbyhistNum, CD8countsbyhistNum ,CD4countsbyhistNum ,FOXP3countsbyhistNum,HLA_ABCPercentagescoresbyhistNum, HLA_ABCHscoresbyhistNum,HLA_C2PercentagescoresbyhistNum, HLA_C2HscoresbyhistNum, HLA_EPercentagescoresbyhistNum, HLA_EHscoresbyhistNum,TAP1PercentagescoresbyhistNum , TAP1HscoresbyhistNum, CSDE1scoresbyhistNum, CSDE1PercentagescoresbyhistNum)
Allmarkermerged <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list)
# remove the orient core
Allmarkermerged <- Allmarkermerged[!Allmarkermerged$HistopathologyNumber == "Orient",]
rownames(Allmarkermerged) <- Allmarkermerged[,1]
#Allmarkermerged <- Allmarkermerged[,-1]
Allmarkermerged <- mutate_all(Allmarkermerged, function(x) as.numeric(as.character(x)))
Allmarkermerged$HistopathologyNumber <- rownames(Allmarkermerged)
omitNAAllmarkerMerged <- na.omit(Allmarkermerged)
write.csv(Allmarkermerged, file = "C:/Users/wp1g19/Documents/TMA analysis/Heatmapper/Quantile cut/TMAanalysis Quantile.csv")
write.csv(omitNAAllmarkerMerged, file = "C:/Users/wp1g19/Documents/TMA analysis/Heatmapper/Quantile cut/NAomitTMAanalysis for heatmapper Quantile.csv")
Allmarkermerged <- Allmarkermerged[,-1]
omitNAAllmarkerMerged <- omitNAAllmarkerMerged[,-1]
#Correlate the counts/scores
setwd("C:/Users/wp1g19/Documents/TMA analysis/R full analysis wd/images for thesis/Quantile_cut")
library("PerformanceAnalytics")
chart.Correlation(Allmarkermerged, histogram=TRUE, pch= "+", method = "pearson")
dev.copy(png,filename="CorrelationQuantileCut.png");
while (dev.cur()>1) dev.off()
chart.Correlation(omitNAAllmarkerMerged, histogram=TRUE, pch= "+", method = "pearson")
dev.copy(png,filename="CorrelationOmitNAQuantileCut.png");
while (dev.cur()>1) dev.off()
#Correaltion heatmap

cormat <- round(cor(as.matrix(omitNAAllmarkerMerged)),2)
library(reshape2)
melted_cormat <- melt(cormat)
head(melted_cormat)
library(ggplot2)
ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()
# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
upper_tri <- get_upper_tri(cormat)
upper_tri
# Melt the correlation matrix
library(reshape2)
melted_cormat <- melt(upper_tri, na.rm = TRUE)

# Heatmap
library(ggplot2)
ggheatmap <- ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson/nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()
#Add correlation coefficients on the heatmap
library(ggplot2)

cor <- Hmisc::rcorr(Allmarkermerged %>% as.matrix(), type = "pearson")
nm = rownames(cor$r)
m = t(combn(nm, 2))
cor$P[5, 13] = 0.049
d = cbind(data.frame(m), R = cor$r[m], P = cor$P[m])
d$label = round(d$R, 2)
d$label[d$P < 0.05] = paste0(d$label[d$P < 0.05], "\n*")
d$X1 = factor(d$X1, nm)
d$X2 = factor(d$X2, rev(nm))

graphics.off()
ggplot(d, aes(X2, X1, fill = R, label = label)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "blue", 
    mid = "white", 
    high = "red", 
    midpoint = .02
  ) +
  geom_text(color = ifelse(d$R > 0.35, "black", "black"),size = 4, lineheight = .8) + 
  theme_bw() +
  coord_equal()+ theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + ggtitle("Correlation of APM stain scores and immune cell density (n = 175)")+ theme(text = element_text(size = 15)) 

pvalues <- as.data.frame(cor$P)  
dev.copy(png,filename="CorrelationOmitNAHeatmapQuantileCut.png");
while (dev.cur()>1) dev.off()

#P adjusted heatmap correlation

library(ggplot2)
library(RcmdrMisc)
library(Hmisc)
cor <- rcorr.adjust(omitNAAllmarkerMerged %>% as.matrix(), type = "pearson")
nm = rownames(cor[["R"]][["r"]])
m = t(combn(nm, 2))
d = cbind(data.frame(m), R = cor[["R"]][["r"]][m], P = cor[["R"]][["P"]][m])
d$label = round(d$R, 2)
d$label[d$P < 0.05] = paste0(d$label[d$P < 0.05], "\n*")
d$X1 = factor(d$X1, nm)
d$X2 = factor(d$X2, rev(nm))

graphics.off()
ggplot(d, aes(X2, X1, fill = R, label = label)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "blue", 
    mid = "white", 
    high = "red", 
    midpoint = .02
  ) +
  geom_text(color = ifelse(d$R > 0.35, "black", "black"),size = 4, lineheight = .8) + 
  theme_bw() +
  coord_equal()+ theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + ggtitle("Correlation of MHC genes and immune cell density (n = 87)")+ theme(text = element_text(size = 15)) 

pvalues <- as.data.frame(cor[["R"]][["P"]])       







regcut <- as.data.frame(omitNAAllmarkerMerged)
regcut1 <- regcut[regcut$HLA_ABCHscore >400,]
regcut2 <- regcut1[regcut1$CSDE1Hscore <20,]
regcut3 <- regcut2[regcut2$`HLA-A` <2000,]
regcut4 <- regcut3[regcut3$`HLA-A` >1,]
regcut5 <- regcut4[regcut4$`HLA-B` <4000,]
regcut6 <- regcut5[regcut5$`HLA-C` <4000,]
regcut$CD8countpermm2
library(ggpubr)
ggplot(regcut, aes(x=CSDE1Hscore, y=HLA_ABCHscore)) + 
  geom_point() +stat_smooth(method = "lm", col = "red") +
  stat_cor(label.x = 3, label.y = 1500) 

#align Clinical data to the counts/scores
rownames(`TMA datasheet final referenced.csv`) <- `TMA datasheet final referenced.csv`[,2]
`TMA datasheet final referenced.csv` <- as.data.frame(`TMA datasheet final referenced.csv`)
Allmarkermerged$Histopathology.Reference.Number <- rownames(Allmarkermerged)
#remove whitespaces from histopathology refertence number
Allmarkermerged %>% 
  mutate(across(where(is.character), str_remove_all, pattern = fixed(" ")))








#save the merged counts/scores
write.csv(Allmarkermerged, file = "C:/Users/wp1g19/Documents/TMA analysis/R full analysis wd/AllmarkersmergedQuantileCut.csv")

#readback the csvs
Allmarkermerged <- read.csv(file = "C:/Users/wp1g19/Documents/TMA analysis/R full analysis wd/AllmarkersmergedQuantileCut.csv", row.names = 1)
Allmarkermerged %>% 
  mutate(across(where(is.character), str_remove_all, pattern = fixed(" ")))
TMACLINmerge <- merge(Allmarkermerged, `TMA datasheet final referenced.csv`, by = "Histopathology.Reference.Number")
#convert survival from days to months
TMACLINmerge <- TMACLINmerge %>% 
  mutate(OS.months = round(OS.survival.time/30.417, digit=0))

TMACLINmerge <- TMACLINmerge %>% 
  mutate(Disease.free.survival.months = round(Disease.free.survival.time/30.417, digit=0))
write.csv(TMACLINmerge, file = "C:/Users/wp1g19/Documents/TMA analysis/R full analysis wd/images for thesis/Quantile_cut/OSkmplotter.csv")


#Set WD to image file


#Overall survival with optimal cutoffs
setwd("C:/Users/wp1g19/Documents/TMA analysis/R full analysis wd/Maximal concordant survival analysis/OS")
# make a cutoff of 120 months
library(survminer)
#adjust the surv_cutpoint formula
surv_cutpoint <- function(data, time = "time", event = "event", variables,
                          minprop = 0.1, maxprop = 0.9, progressbar = TRUE)
{
  if(!inherits(data, "data.frame"))
    stop("data should be an object of class data.frame")
  data <- as.data.frame(data)
  if(!all(c(time, event) %in% colnames(data)))
    stop("Specify correct column names containing time and event values.")
  if(!all(variables %in% colnames(data)))
    stop("Some variables are not found in the data: ",
         paste(setdiff(variables, colnames(data)), collapse =", "))
  
  not_numeric_vars <- .get_not_numeric_vars(data[, variables, drop = FALSE])
  variables <- setdiff(variables, not_numeric_vars) # keep only numeric variables
  if(length(variables)==0) stop("At least, one numeric variables required.")
  
  nvar <- length(variables)
  if(nvar <= 5) progressbar <- FALSE
  if(progressbar) pb <- utils::txtProgressBar(min = 0, max = nvar, style = 3)
  surv_data <- data.frame(time = data[, time], event = data[, event])
  res <- list()
  for (i in 1:nvar){
    var_i <- variables[i]
    surv_data$var <- data[, var_i]
    max_stat_i <- maxstat::maxstat.test(survival::Surv(time, event) ~ var, data = surv_data,
                                        smethod = "LogRank", pmethod="none",
                                        minprop = minprop, maxprop = maxprop,
                                        alpha = alpha)
    res[[var_i]] <- max_stat_i
    if(progressbar) utils::setTxtProgressBar(pb, i)
  }
  colnames(surv_data) <- c(time, event)
  res$data <- cbind.data.frame(surv_data[, 1:2, drop = FALSE], data[, variables, drop = FALSE])
  res$minprop <- minprop
  if(!is.null(not_numeric_vars)) res$not_numeric <- data[, not_numeric_vars, drop = FALSE]
  res <- structure(res, class = c("list", "surv_cutpoint"))
  res$cutpoint <- summary(res)
  res
}

.get_not_numeric_vars <- function(data_frame){
  is_numeric <- sapply(data_frame, is.numeric)
  if(sum(!is_numeric) == 0) res = NULL
  else res <- colnames(data_frame[, !is_numeric, drop = FALSE])
  res
}





################# CD3 survival


#make the optimal cut
CD3countpermm2.surv_count.cut <- surv_cutpoint(
  TMACLINmerge,
  time = "OS.months",
  event = "Survival.status",
  variables = "CD3countpermm2", minprop = 0.45, maxprop = 0.65
)

cd3cutsum <- as.data.frame(summary(CD3countpermm2.surv_count.cut))

#Plot the cutpoint 

# palette = "npg" (nature publishing group), see ?ggpubr::ggpar
plot(CD3countpermm2.surv_count.cut, "CD3countpermm2", palette = "npg", bins = 152, xlim=c(400, 700))
dev.copy(png,filename="CD3.surv_count.cut.png");
while (dev.cur()>1) dev.off()

#categorise the cutpoint 
CD3.surv_count.cat <- surv_categorize(CD3countpermm2.surv_count.cut) 

library(survival)
library(RTCGA)
fit <- survfit(Surv(OS.months, Survival.status) ~ CD3countpermm2,
               data = CD3.surv_count.cat)
# extract p value
cd3pvalue <- surv_pvalue(fit)$pval.txt
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,160), 
  xlab = ("Months"),# present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 20,    # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = TRUE # show bars instead of names in text annotations
  # in legend of risk table
)  
#save the image file
dev.copy(png,filename="CD3OSKm.png");
while (dev.cur()>1) dev.off()


########### CD4 survival


#make the optimal cut
CD4countpermm2.surv_count.cut <- surv_cutpoint(
  TMACLINmerge,
  time = "OS.months",
  event = "Survival.status",
  variables = "CD4countpermm2", minprop = 0.3, maxprop = 0.4
)

cd4cutsum <- as.data.frame(summary(CD4countpermm2.surv_count.cut))

#Plot the cutpoint 

# palette = "npg" (nature publishing group), see ?ggpubr::ggpar
plot(CD4countpermm2.surv_count.cut, "CD4countpermm2", palette = "npg", bins = 166)
dev.copy(png,filename="CD4.surv_count.cut.png");
while (dev.cur()>1) dev.off()

#categorise the cutpoint 
CD4.surv_count.cat <- surv_categorize(CD4countpermm2.surv_count.cut) 

library(survival)
library(RTCGA)
fit <- survfit(Surv(OS.months, Survival.status) ~ CD4countpermm2,
               data = CD4.surv_count.cat)
# extract p value
cd4pvalue <- surv_pvalue(fit)$pval.txt
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,180), 
  xlab = ("Months"),# present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 20,    # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = TRUE # show bars instead of names in text annotations
  # in legend of risk table
)  
#save the image file
dev.copy(png,filename="CD4OSKm.png");
while (dev.cur()>1) dev.off()

########### CD8 survival


#make the optimal cut
CD8countpermm2.surv_count.cut <- surv_cutpoint(
  TMACLINmerge,
  time = "OS.months",
  event = "Survival.status",
  variables = "CD8countpermm2", minprop = 0.2, maxprop = 0.8
)

cd8cutsum <- as.data.frame(summary(CD8countpermm2.surv_count.cut))

#Plot the cutpoint 

# palette = "npg" (nature publishing group), see ?ggpubr::ggpar
plot(CD8countpermm2.surv_count.cut, "CD8countpermm2", palette = "npg", bins = 166)
dev.copy(png,filename="CD8.surv_count.cut.png");
while (dev.cur()>1) dev.off()

#categorise the cutpoint 
CD8.surv_count.cat <- surv_categorize(CD8countpermm2.surv_count.cut) 

library(survival)
library(RTCGA)
fit <- survfit(Surv(OS.months, Survival.status) ~ CD8countpermm2,
               data = CD8.surv_count.cat)
# extract p value
cd8pvalue <- surv_pvalue(fit)$pval.txt
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,180), 
  xlab = ("Months"),# present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 20,    # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = TRUE # show bars instead of names in text annotations
  # in legend of risk table
)  
#save the image file
dev.copy(png,filename="CD8OSKm.png");
while (dev.cur()>1) dev.off()

########### FOXP3 survival


#make the optimal cut
FOXP3countpermm2.surv_count.cut <- surv_cutpoint(
  TMACLINmerge,
  time = "OS.months",
  event = "Survival.status",
  variables = "FOXP3countpermm2", minprop = 0.2, maxprop = 0.8
)

foxp3cutsum <- as.data.frame(summary(FOXP3countpermm2.surv_count.cut))

#Plot the cutpoint 

# palette = "npg" (nature publishing group), see ?ggpubr::ggpar
plot(FOXP3countpermm2.surv_count.cut, "FOXP3countpermm2", palette = "npg", bins = 166)
dev.copy(png,filename="FOXP3.surv_count.cut.png");
while (dev.cur()>1) dev.off()

#categorise the cutpoint 
FOXP3.surv_count.cat <- surv_categorize(FOXP3countpermm2.surv_count.cut) 

library(survival)
library(RTCGA)
fit <- survfit(Surv(OS.months, Survival.status) ~ FOXP3countpermm2,
               data = FOXP3.surv_count.cat)
# extract p value
foxp3pvalue <- surv_pvalue(fit)$pval.txt
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,180), 
  xlab = ("Months"),# present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 20,    # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = TRUE # show bars instead of names in text annotations
  # in legend of risk table
)  
#save the image file
dev.copy(png,filename="FOXP3OSKm.png");
while (dev.cur()>1) dev.off()

########### HLA-ABC percent survival

#make the optimal cut
HLA_ABCpercentagescore.surv_count.cut <- surv_cutpoint(
  TMACLINmerge,
  time = "OS.months",
  event = "Survival.status",
  variables = "HLA_ABCpercentagescore", minprop = 0.1, maxprop = 0.8
)

HLA_ABC_perccutsum <- as.data.frame(summary(HLA_ABCpercentagescore.surv_count.cut))

#Plot the cutpoint 

# palette = "npg" (nature publishing group), see ?ggpubr::ggpar
plot(HLA_ABCpercentagescore.surv_count.cut, "HLA_ABCpercentagescore", palette = "npg", bins = 166)
dev.copy(png,filename="HLA_ABCpercentagescore.surv_count.cut.png");
while (dev.cur()>1) dev.off()

#categorise the cutpoint 
HLA_ABCpercentagescore.surv_count.cat <- surv_categorize(HLA_ABCpercentagescore.surv_count.cut) 

library(survival)
library(RTCGA)
fit <- survfit(Surv(OS.months, Survival.status) ~ HLA_ABCpercentagescore,
               data = HLA_ABCpercentagescore.surv_count.cat)
# extract p value
HLAABCperpvalue <- surv_pvalue(fit)$pval.txt
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,60), 
  xlab = ("Months"),# present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 20,    # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(),  font.x = 25, font.y = 25, font.legend = 25, font.caption = 25,  font.tickslab = 12,  # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = TRUE, surv.median.line = "hv" # show bars instead of names in text annotations
  # in legend of risk table
)  + 
  guides(colour = guide_legend(nrow = 2))
#save the image file
dev.copy(png,filename="HLA_ABCpercentagescoreOSKm.png");
while (dev.cur()>1) dev.off()


############# HLA_ABC H score survival
#make the optimal cut
HLA_ABCHscore.surv_count.cut <- surv_cutpoint(
  TMACLINmerge,
  time = "OS.months",
  event = "Survival.status",
  variables = "HLA_ABCHscore", minprop = 0.2, maxprop = 0.4
)

HLA_ABC_HScutsum <- as.data.frame(summary(HLA_ABCHscore.surv_count.cut))

#Plot the cutpoint 

# palette = "npg" (nature publishing group), see ?ggpubr::ggpar
plot(HLA_ABCHscore.surv_count.cut, "HLA_ABCHscore", palette = "npg", bins = 166)
dev.copy(png,filename="HLA_ABCHscore.surv_count.cut.png");
while (dev.cur()>1) dev.off()

#categorise the cutpoint 
HLA_ABCHscore.surv_count.cat <- surv_categorize(HLA_ABCHscore.surv_count.cut) 

library(survival)
library(RTCGA)
fit <- survfit(Surv(OS.months, Survival.status) ~ HLA_ABCHscore,
               data = HLA_ABCHscore.surv_count.cat)
# extract p value
HLAABCHSpvalue <- surv_pvalue(fit)$pval.txt
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  risk.table = FALSE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,60), 
  xlab = ("Months"),# present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 20,    # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(),   font.x = 25, font.y = 25, font.legend = 25, font.caption = 25,  font.tickslab = 20, pval.size = 10, # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = TRUE, surv.median.line = "hv" # show bars instead of names in text annotations
  # in legend of risk table
)   + 
  guides(colour = guide_legend(nrow = 2))
#save the image file
dev.copy(png,filename="HLA_ABCHscoreOSKm.png");
while (dev.cur()>1) dev.off()



############# HLA_ABC H score survival
#make the optimal cut
HLA_C2percentagescore.surv_count.cut <- surv_cutpoint(
  TMACLINmerge,
  time = "OS.months",
  event = "Survival.status",
  variables = "HLA_C2positive.", minprop = 0.2, maxprop = 0.8
)

HLA_C2_perccutsum <- as.data.frame(summary(HLA_C2percentagescore.surv_count.cut))

#Plot the cutpoint 

# palette = "npg" (nature publishing group), see ?ggpubr::ggpar
plot(HLA_C2percentagescore.surv_count.cut, "HLA_C2positive.", palette = "npg", bins = 166)
dev.copy(png,filename="HLA_C2percentagescore.surv_count.cut.png");
while (dev.cur()>1) dev.off()

#categorise the cutpoint 
HLA_C2percentagescore.surv_count.cat <- surv_categorize(HLA_C2percentagescore.surv_count.cut) 

library(survival)
library(RTCGA)
fit <- survfit(Surv(OS.months, Survival.status) ~ HLA_C2positive.,
               data = HLA_C2percentagescore.surv_count.cat)
# extract p value
HLAC2perpvalue <- surv_pvalue(fit)$pval.txt
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,180), 
  xlab = ("Months"),# present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 20,    # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = TRUE # show bars instead of names in text annotations
  # in legend of risk table
)  
#save the image file
dev.copy(png,filename="HLA_C2percentagescoreOSKm.png");
while (dev.cur()>1) dev.off()

############# HLA_C2 H score survival
#make the optimal cut
HLA_C2Hscore.surv_count.cut <- surv_cutpoint(
  TMACLINmerge,
  time = "OS.months",
  event = "Survival.status",
  variables = "HLA_C2Hscore", minprop = 0.2, maxprop = 0.8
)

HLA_C2_HScutsum <- as.data.frame(summary(HLA_C2Hscore.surv_count.cut))

#Plot the cutpoint 

# palette = "npg" (nature publishing group), see ?ggpubr::ggpar
plot(HLA_C2Hscore.surv_count.cut, "HLA_C2Hscore", palette = "npg", bins = 166)
dev.copy(png,filename="HLA_C2Hscore.surv_count.cut.png");
while (dev.cur()>1) dev.off()

#categorise the cutpoint 
HLA_C2Hscore.surv_count.cat <- surv_categorize(HLA_C2Hscore.surv_count.cut) 

library(survival)
library(RTCGA)
fit <- survfit(Surv(OS.months, Survival.status) ~ HLA_C2Hscore,
               data = HLA_C2Hscore.surv_count.cat)
# extract p value
HLAC2HSpvalue <- surv_pvalue(fit)$pval.txt
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  risk.table = FALSE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,60), 
  xlab = ("Months"),# present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 20,    # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(),font.x = 25, font.y = 25, font.legend = 25, font.caption = 25,  font.tickslab = 20, pval.size = 10, # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = TRUE, surv.median.line = "hv" # show bars instead of names in text annotations
  # in legend of risk table
)  + 
  guides(colour = guide_legend(nrow = 2))
#save the image file
dev.copy(png,filename="HLA_C2HscoreOSKm.png");
while (dev.cur()>1) dev.off()

############# HLA_Epercentagescore survival
#make the optimal cut
HLA_Epercentagescore.surv_count.cut <- surv_cutpoint(
  TMACLINmerge,
  time = "OS.months",
  event = "Survival.status",
  variables = "HLA_Epositive.", minprop = 0.45,maxprop = 0.55
)

HLA_E_perccutsum <- as.data.frame(summary(HLA_Epercentagescore.surv_count.cut))

#Plot the cutpoint 

# palette = "npg" (nature publishing group), see ?ggpubr::ggpar
plot(HLA_Epercentagescore.surv_count.cut, "HLA_Epositive.", palette = "npg", bins = 166)
dev.copy(png,filename="HLA_Epercentagescore.surv_count.cut.png");
while (dev.cur()>1) dev.off()

#categorise the cutpoint 
HLA_Epercentagescore.surv_count.cat <- surv_categorize(HLA_Epercentagescore.surv_count.cut) 

library(survival)
library(RTCGA)
fit <- survfit(Surv(OS.months, Survival.status) ~ HLA_Epositive.,
               data = HLA_Epercentagescore.surv_count.cat)
# extract p value
HLAEperpvalue <- surv_pvalue(fit)$pval.txt
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,180), 
  xlab = ("Months"),# present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 20,    # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = TRUE # show bars instead of names in text annotations
  # in legend of risk table
)  
#save the image file
dev.copy(png,filename="HLA_EpercentagescoreOSKm.png");
while (dev.cur()>1) dev.off()

############# HLA_EHscore survival
#make the optimal cut
HLA_EHscore.surv_count.cut <- surv_cutpoint(
  TMACLINmerge,
  time = "OS.months",
  event = "Survival.status",
  variables = "HLA_EHscore"
)

HLA_E_HScutsum <- as.data.frame(summary(HLA_EHscore.surv_count.cut))

#Plot the cutpoint 

# palette = "npg" (nature publishing group), see ?ggpubr::ggpar
plot(HLA_EHscore.surv_count.cut, "HLA_EHscore", palette = "npg", bins = 166)
dev.copy(png,filename="HLA_EHscore.surv_count.cut.png");
while (dev.cur()>1) dev.off()

#categorise the cutpoint 
HLA_EHscore.surv_count.cat <- surv_categorize(HLA_EHscore.surv_count.cut) 

library(survival)
library(RTCGA)
fit <- survfit(Surv(OS.months, Survival.status) ~ HLA_EHscore,
               data = HLA_EHscore.surv_count.cat)
# extract p value
HLAEHSpvalue <- surv_pvalue(fit)$pval.txt
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  risk.table = FALSE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,60), 
  xlab = ("Months"),# present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 20,    # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(),  font.x = 25, font.y = 25, font.legend = 25, font.caption = 25,  font.tickslab = 20, pval.size = 10, # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = TRUE, surv.median.line = "hv" # show bars instead of names in text annotations
  # in legend of risk table
)  + 
  guides(colour = guide_legend(nrow = 2))
#save the image file
dev.copy(png,filename="HLA_EHscoreOSKm.png");
while (dev.cur()>1) dev.off()

############# TAP1percentagescore survival
#make the optimal cut
TAP1percentagescore.surv_count.cut <- surv_cutpoint(
  TMACLINmerge,
  time = "OS.months",
  event = "Survival.status",
  variables = "TAP1positive.", minprop = 0.2, maxprop = 0.7
)

TAP1_perccutsum <- as.data.frame(summary(TAP1percentagescore.surv_count.cut))

#Plot the cutpoint 

# palette = "npg" (nature publishing group), see ?ggpubr::ggpar
plot(TAP1percentagescore.surv_count.cut, "TAP1positive.", palette = "npg", bins = 166)
 dev.copy(png,filename="TAP1percentagescore.surv_count.cut.png");
while (dev.cur()>1) dev.off()

#categorise the cutpoint 
TAP1percentagescore.surv_count.cat <- surv_categorize(TAP1percentagescore.surv_count.cut) 

library(survival)
library(RTCGA)
fit <- survfit(Surv(OS.months, Survival.status) ~ TAP1positive.,
               data = TAP1percentagescore.surv_count.cat)
# extract p value
TAP1perpvalue <- surv_pvalue(fit)$pval.txt
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,180), 
  xlab = ("Months"),# present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 20,    # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = TRUE # show bars instead of names in text annotations
  # in legend of risk table
)  
#save the image file
dev.copy(png,filename="TAP1percentagescoreOSKm.png");
while (dev.cur()>1) dev.off()

############# TAP1Hscore survival
#make the optimal cut
TAP1Hscore.surv_count.cut <- surv_cutpoint(
  TMACLINmerge,
  time = "OS.months",
  event = "Survival.status",
  variables = "TAP1Hscore", minprop = 0.25, maxprop = 0.75
)

TAP1_HScutsum <- as.data.frame(summary(TAP1Hscore.surv_count.cut))

#Plot the cutpoint 

# palette = "npg" (nature publishing group), see ?ggpubr::ggpar
plot(TAP1Hscore.surv_count.cut, "TAP1Hscore", palette = "npg", bins = 166)
dev.copy(png,filename="TAP1Hscore.surv_count.cut.png");
while (dev.cur()>1) dev.off()

#categorise the cutpoint 
TAP1Hscore.surv_count.cat <- surv_categorize(TAP1Hscore.surv_count.cut) 

library(survival)
library(RTCGA)
fit <- survfit(Surv(OS.months, Survival.status) ~ TAP1Hscore,
               data = TAP1Hscore.surv_count.cat)
# extract p value
TAP1HSpvalue <- surv_pvalue(fit)$pval.txt

ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,180), 
  xlab = ("Months"),# present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 20,    # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = TRUE # show bars instead of names in text annotations
  # in legend of risk table
)  
#save the image file
dev.copy(png,filename="TAP1HscoreOSKm.png");
while (dev.cur()>1) dev.off()


#bind the cutpoint point summaries
cutpointsums <- rbind(cd3cutsum,
      cd4cutsum,
      cd8cutsum,
      foxp3cutsum,
      HLA_ABC_perccutsum,
      HLA_ABC_HScutsum,
      HLA_C2_perccutsum,
      HLA_C2_HScutsum,
      HLA_E_perccutsum,
      HLA_E_HScutsum,
      TAP1_perccutsum,
      TAP1_HScutsum)
pvalues <- as.data.frame(rbind(cd3pvalue,
                               cd4pvalue,
                               cd8pvalue,
                               foxp3pvalue,
                               HLAABCperpvalue,
                               HLAABCHSpvalue,
                               HLAC2perpvalue,
                               HLAC2HSpvalue,
                               HLAEperpvalue,
                               HLAEHSpvalue,
                               TAP1perpvalue,
                               TAP1HSpvalue))
cutpointpvaluesums <- as.data.frame(cutpointsums)
cutpointpvaluesums$pvalueOS <- pvalues
write.csv(cutpointpvaluesums, file = "C:/Users/wp1g19/Documents/TMA analysis/R full analysis wd/images for thesis/Quantile_cut/quantcutpointpvalues.csv")

#Univariate CoxPH
# Place the optimal cutpoint catagories into a single dataframe

#make the optimal cut on dataframe with NAs removed
TMACLINmergeCD3 <- TMACLINmerge[,c("Histopathology.Reference.Number","OS.months", "Survival.status", "CD3countpermm2")]
CD3countpermm2.surv_count.cut.histo <- surv_cutpoint(
  TMACLINmergeCD3,
  time = "OS.months",
  event = "Survival.status",
  variables = "CD3countpermm2", minprop = 0.45, maxprop = 0.65
)

summary(CD3countpermm2.surv_count.cut)


#categorise the cutpoint 
CD3.surv_count.cat.histo <- surv_categorize(CD3countpermm2.surv_count.cut.histo) 
CD3.surv_count.cat.histo$Histopathology.Reference.Number <- TMACLINmergeCD3$Histopathology.Reference.Number


TMACLINmergeCD4 <- TMACLINmerge[,c("Histopathology.Reference.Number","OS.months", "Survival.status", "CD4countpermm2")]
CD4countpermm2.surv_count.cut.histo <- surv_cutpoint(
  TMACLINmergeCD4,
  time = "OS.months",
  event = "Survival.status",
  variables = "CD4countpermm2", minprop = 0.3, maxprop = 0.4
)

summary(CD4countpermm2.surv_count.cut.histo)


#categorise the cutpoint 
CD4.surv_count.cat.histo <- surv_categorize(CD4countpermm2.surv_count.cut.histo) 
CD4.surv_count.cat.histo$Histopathology.Reference.Number <- TMACLINmergeCD4$Histopathology.Reference.Number



TMACLINmergeCD8 <- TMACLINmerge[,c("Histopathology.Reference.Number","OS.months", "Survival.status", "CD8countpermm2")]
CD8countpermm2.surv_count.cut.histo <- surv_cutpoint(
  TMACLINmergeCD8,
  time = "OS.months",
  event = "Survival.status",
  variables = "CD8countpermm2", minprop = 0.1, maxprop = 0.7
)

summary(CD8countpermm2.surv_count.cut.histo)


#categorise the cutpoint 
CD8.surv_count.cat.histo <- surv_categorize(CD8countpermm2.surv_count.cut.histo) 
CD8.surv_count.cat.histo$Histopathology.Reference.Number <- TMACLINmergeCD8$Histopathology.Reference.Number


TMACLINmergeFOXP3 <- TMACLINmerge[,c("Histopathology.Reference.Number","OS.months", "Survival.status", "FOXP3countpermm2")]
FOXP3countpermm2.surv_count.cut.histo <- surv_cutpoint(
  TMACLINmergeFOXP3,
  time = "OS.months",
  event = "Survival.status",
  variables = "FOXP3countpermm2", minprop = 0.2, maxprop = 0.8
)

summary(FOXP3countpermm2.surv_count.cut.histo)


#categorise the cutpoint 
FOXP3.surv_count.cat.histo <- surv_categorize(FOXP3countpermm2.surv_count.cut.histo) 
FOXP3.surv_count.cat.histo$Histopathology.Reference.Number <- TMACLINmergeFOXP3$Histopathology.Reference.Number



TMACLINmergeHLA_ABCpercentagescore <- TMACLINmerge[,c("Histopathology.Reference.Number","OS.months", "Survival.status", "HLA_ABCpercentagescore")]
HLA_ABCpercentagescore.surv_count.cut.histo <- surv_cutpoint(
  TMACLINmergeHLA_ABCpercentagescore,
  time = "OS.months",
  event = "Survival.status",
  variables = "HLA_ABCpercentagescore", minprop = 0.2, maxprop = 0.4
)

summary(HLA_ABCpercentagescore.surv_count.cut.histo)


#categorise the cutpoint 
HLA_ABCpercentagescore.surv_count.cat.histo <- surv_categorize(HLA_ABCpercentagescore.surv_count.cut.histo) 
HLA_ABCpercentagescore.surv_count.cat.histo$Histopathology.Reference.Number <- TMACLINmergeHLA_ABCpercentagescore$Histopathology.Reference.Number



TMACLINmergeHLA_ABCHscore <- TMACLINmerge[,c("Histopathology.Reference.Number","OS.months", "Survival.status", "HLA_ABCHscore")]
HLA_ABCHscore.surv_count.cut.histo <- surv_cutpoint(
  TMACLINmergeHLA_ABCHscore,
  time = "OS.months",
  event = "Survival.status",
  variables = "HLA_ABCHscore", minprop = 0.2, maxprop = 0.4
)

summary(HLA_ABCHscore.surv_count.cut.histo)


#categorise the cutpoint 
HLA_ABCHscore.surv_count.cat.histo <- surv_categorize(HLA_ABCHscore.surv_count.cut.histo) 
HLA_ABCHscore.surv_count.cat.histo$Histopathology.Reference.Number <- TMACLINmergeHLA_ABCHscore$Histopathology.Reference.Number

TMACLINmergeHLA_C2percentagescore <- TMACLINmerge[,c("Histopathology.Reference.Number","OS.months", "Survival.status", "HLA_C2positive.")]
HLA_C2percentagescore.surv_count.cut.histo <- surv_cutpoint(
  TMACLINmergeHLA_C2percentagescore,
  time = "OS.months",
  event = "Survival.status",
  variables = "HLA_C2positive.", minprop = 0.2, maxprop = 0.8
)

summary(HLA_C2percentagescore.surv_count.cut.histo)


#categorise the cutpoint 
HLA_C2percentagescore.surv_count.cat.histo <- surv_categorize(HLA_C2percentagescore.surv_count.cut.histo) 
HLA_C2percentagescore.surv_count.cat.histo$Histopathology.Reference.Number <- TMACLINmergeHLA_C2percentagescore$Histopathology.Reference.Number


TMACLINmergeHLA_C2Hscore <- TMACLINmerge[,c("Histopathology.Reference.Number","OS.months", "Survival.status", "HLA_C2Hscore")]
HLA_C2Hscore.surv_count.cut.histo <- surv_cutpoint(
  TMACLINmergeHLA_C2Hscore,
  time = "OS.months",
  event = "Survival.status",
  variables = "HLA_C2Hscore", minprop = 0.2, maxprop = 0.8
)

summary(HLA_C2Hscore.surv_count.cut.histo)


#categorise the cutpoint 
HLA_C2Hscore.surv_count.cat.histo <- surv_categorize(HLA_C2Hscore.surv_count.cut.histo) 
HLA_C2Hscore.surv_count.cat.histo$Histopathology.Reference.Number <- TMACLINmergeHLA_C2Hscore$Histopathology.Reference.Number


TMACLINmergeHLA_Epercentagescore <- TMACLINmerge[,c("Histopathology.Reference.Number","OS.months", "Survival.status", "HLA_Epositive.")]
HLA_Epercentagescore.surv_count.cut.histo <- surv_cutpoint(
  TMACLINmergeHLA_Epercentagescore,
  time = "OS.months",
  event = "Survival.status",
  variables = "HLA_Epositive.", minprop = 0.45, maxprop = 0.55
)

summary(HLA_Epercentagescore.surv_count.cut.histo)


#categorise the cutpoint 
HLA_Epercentagescore.surv_count.cat.histo <- surv_categorize(HLA_Epercentagescore.surv_count.cut.histo) 
HLA_Epercentagescore.surv_count.cat.histo$Histopathology.Reference.Number <- TMACLINmergeHLA_Epercentagescore$Histopathology.Reference.Number

TMACLINmergeHLA_EHscore <- TMACLINmerge[,c("Histopathology.Reference.Number","OS.months", "Survival.status", "HLA_EHscore")]
HLA_EHscore.surv_count.cut.histo <- surv_cutpoint(
  TMACLINmergeHLA_EHscore,
  time = "OS.months",
  event = "Survival.status",
  variables = "HLA_EHscore"
)

summary(HLA_EHscore.surv_count.cut.histo)


#categorise the cutpoint 
HLA_EHscore.surv_count.cat.histo <- surv_categorize(HLA_EHscore.surv_count.cut.histo) 
HLA_EHscore.surv_count.cat.histo$Histopathology.Reference.Number <- TMACLINmergeHLA_EHscore$Histopathology.Reference.Number


TMACLINmergeTAP1Hscore <- TMACLINmerge[,c("Histopathology.Reference.Number","OS.months", "Survival.status", "TAP1Hscore")]
TAP1Hscore.surv_count.cut.histo <- surv_cutpoint(
  TMACLINmergeTAP1Hscore,
  time = "OS.months",
  event = "Survival.status",
  variables = "TAP1Hscore", minprop = 0.25, maxprop = 0.75
)

summary(TAP1Hscore.surv_count.cut.histo)


#categorise the cutpoint 
TAP1Hscore.surv_count.cat.histo <- surv_categorize(TAP1Hscore.surv_count.cut.histo) 
TAP1Hscore.surv_count.cat.histo$Histopathology.Reference.Number <- TMACLINmergeTAP1Hscore$Histopathology.Reference.Number

TMACLINmergeTAP1percentagescore <- TMACLINmerge[,c("Histopathology.Reference.Number","OS.months", "Survival.status", "TAP1positive.")]
TAP1percentagescore.surv_count.cut.histo <- surv_cutpoint(
  TMACLINmergeTAP1percentagescore,
  time = "OS.months",
  event = "Survival.status",
  variables = "TAP1positive.", minprop = 0.2, maxprop = 0.7
)

summary(TAP1percentagescore.surv_count.cut.histo)


#categorise the cutpoint 
TAP1percentagescore.surv_count.cat.histo <- surv_categorize(TAP1percentagescore.surv_count.cut.histo) 
TAP1percentagescore.surv_count.cat.histo$Histopathology.Reference.Number <- TMACLINmergeTAP1percentagescore$Histopathology.Reference.Number

TMACLINmergeCSDE1percentagescore <- TMACLINmerge[,c("Histopathology.Reference.Number","OS.months", "Survival.status", "CSDE1positive.")]
CSDE1percentagescore.surv_count.cut.histo <- surv_cutpoint(
  TMACLINmergeCSDE1percentagescore,
  time = "OS.months",
  event = "Survival.status",
  variables = "CSDE1positive.", minprop = 0.2, maxprop = 0.7
)

summary(CSDE1percentagescore.surv_count.cut.histo)


#categorise the cutpoint 
CSDE1percentagescore.surv_count.cat.histo <- surv_categorize(CSDE1percentagescore.surv_count.cut.histo) 
CSDE1percentagescore.surv_count.cat.histo$Histopathology.Reference.Number <- TMACLINmergeCSDE1percentagescore$Histopathology.Reference.Number

TMACLINmergeCSDE1percentagescore <- TMACLINmerge[,c("Histopathology.Reference.Number","OS.months", "Survival.status", "CSDE1positive.")]
CSDE1percentagescore.surv_count.cut.histo <- surv_cutpoint(
  TMACLINmergeCSDE1percentagescore,
  time = "OS.months",
  event = "Survival.status",
  variables = "CSDE1positive.", minprop = 0.2, maxprop = 0.7
)

summary(CSDE1percentagescore.surv_count.cut.histo)


#categorise the cutpoint 
CSDE1percentagescore.surv_count.cat.histo <- surv_categorize(CSDE1percentagescore.surv_count.cut.histo) 
CSDE1percentagescore.surv_count.cat.histo$Histopathology.Reference.Number <- TMACLINmergeCSDE1percentagescore$Histopathology.Reference.Number

TMACLINmergeCSDE1Hscore <- TMACLINmerge[,c("Histopathology.Reference.Number","OS.months", "Survival.status", "CSDE1Hscore")]
CSDE1Hscore.surv_count.cut.histo <- surv_cutpoint(
  TMACLINmergeCSDE1Hscore,
  time = "OS.months",
  event = "Survival.status",
  variables = "CSDE1Hscore", minprop = 0.25, maxprop = 0.75
)

summary(CSDE1Hscore.surv_count.cut.histo)


#categorise the cutpoint 
CSDE1Hscore.surv_count.cat.histo <- surv_categorize(CSDE1Hscore.surv_count.cut.histo) 
CSDE1Hscore.surv_count.cat.histo$Histopathology.Reference.Number <- TMACLINmergeCSDE1Hscore$Histopathology.Reference.Number

#Merge the data
df_list <- list(CD3.surv_count.cat.histo, 
                CD4.surv_count.cat.histo,
                CD8.surv_count.cat.histo,
                FOXP3.surv_count.cat.histo, 
                HLA_ABCpercentagescore.surv_count.cat.histo, 
                HLA_ABCHscore.surv_count.cat.histo, 
                HLA_Epercentagescore.surv_count.cat.histo, 
                HLA_EHscore.surv_count.cat.histo,
                HLA_C2percentagescore.surv_count.cat.histo,
                HLA_C2Hscore.surv_count.cat.histo,
                TAP1percentagescore.surv_count.cat.histo,
                TAP1Hscore.surv_count.cat.histo,
                CSDE1Hscore.surv_count.cat.histo,
                CSDE1percentagescore.surv_count.cat.histo)
COXphoptimalcutoffs <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list)
clincovariates <- c("Histopathology.Reference.Number",
                    "Age",
                    "Sex",
                    "pT",
                    "pN",
                    "pM",
                    "LymphaticInvasion",
                    "VascularInvasion",
                    "Perineural_Invasion",
                    "TumourGradingDifferentiation",
                    "R0.American.",
                    "MandardScore",
                    "AliveorDead")

Multivariateclinicalcovariates <- TMACLINmerge[,clincovariates]


df_list <- list(COXphoptimalcutoffs, Multivariateclinicalcovariates)

COXphoptimalcutoffsClinMerge <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list)
write.csv(COXphoptimalcutoffsClinMerge, file = "C:/Users/wp1g19/Documents/TMA analysis/R full analysis wd/COXphoptimalcutoffsClinMerge.csv")
COXphoptimalcutoffsClinMerge <- read.csv(file = "C:/Users/wp1g19/Documents/TMA analysis/R full analysis wd/COXphoptimalcutoffsClinMerge.csv", header = TRUE)
COXphoptimalcutoffsClinMerge$CSDE1.percentage.positive
covariates <- c("CD3countpermm2",
                "CD4countpermm2",
                "CD8countpermm2", 
                "FOXP3countpermm2", 
                "HLA_ABC.percentage.positive",
                "HLA_ABCHscore",
                "HLA_E.percentage.positive",
                "HLA_EHscore",
                "HLA_C2.percentage.positive",
                "HLA_C2Hscore",
                "TAP1.percentage.positive",
                "TAP1Hscore",
                "CSDE1.percentage.positive",
                "CSDE1Hscore")
                                      
#Forest plot coxPH analysis OS                                  
library(forestmodel)
#relevel for all positive HR
COXphoptimalcutoffsClinMerge$HLA_ABC.percentage.positive <- as.factor(COXphoptimalcutoffsClinMerge$HLA_ABC.percentage.positive)
COXphoptimalcutoffsClinMerge$HLA_ABC.percentage.positive = relevel(COXphoptimalcutoffsClinMerge$HLA_ABC.percentage.positive, ref = "low")
COXphoptimalcutoffsClinMerge$TAP1.percentage.positive <- as.factor(COXphoptimalcutoffsClinMerge$TAP1.percentage.positive)
COXphoptimalcutoffsClinMerge$TAP1.percentage.positive = relevel(COXphoptimalcutoffsClinMerge$TAP1.percentage.positive, ref = "low")
COXphoptimalcutoffsClinMerge$TAP1Hscore <- as.factor(COXphoptimalcutoffsClinMerge$TAP1Hscore)
COXphoptimalcutoffsClinMerge$TAP1Hscore = relevel(COXphoptimalcutoffsClinMerge$TAP1Hscore, ref = "low")
COXphoptimalcutoffsClinMerge$CSDE1.percentage.positive <- as.factor(COXphoptimalcutoffsClinMerge$CSDE1.percentage.positive)
COXphoptimalcutoffsClinMerge$CSDE1.percentage.positive = relevel(COXphoptimalcutoffsClinMerge$CSDE1.percentage.positive, ref = "low")
COXphoptimalcutoffsClinMerge$CSDE1Hscore <- as.factor(COXphoptimalcutoffsClinMerge$CSDE1Hscore)
COXphoptimalcutoffsClinMerge$CSDE1Hscore = relevel(COXphoptimalcutoffsClinMerge$CSDE1Hscore, ref = "low")
COXphoptimalcutoffsClinMerge$HLA_C2Hscore <- as.factor(COXphoptimalcutoffsClinMerge$HLA_C2Hscore)
COXphoptimalcutoffsClinMerge$HLA_C2Hscore = relevel(COXphoptimalcutoffsClinMerge$HLA_C2Hscore, ref = "high")
univ_formulas <- sapply(covariates, function(x) as.formula(paste('Surv(OS.months, Survival.status)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = COXphoptimalcutoffsClinMerge)})


forest_model(model_list = univ_models,
             covariates = covariates,
             merge_models = T,
             recalculate_height = 15,
             recalculate_width = 3)

library(grid)
grid.text("Univariate OS", .3, 0.975, gp=gpar(cex=3, fontsize = 6), check.overlap = T)
dev.copy(png,filename="UnivariateCoxPhForest.png");
while (dev.cur()>1) dev.off()




# Multivariate CoxPH analysis
#merge the optimal coxph cutoff with the full clinical data
library(tidyverse)
library(tidytidbits)
library(survivalAnalysis)

#Create a dataframe with the clinical characeristics to assess in multivariate

clincovariates <- c("Histopathology.Reference.Number",
                    "Age",
                    "Sex",
                    "pT",
                    "pN",
                    "pM",
                    "LymphaticInvasion",
                    "VascularInvasion",
                    "Perineural_Invasion",
                    "TumourGradingDifferentiation",
                    "R0.American.",
                    "MandardScore",
                    "AliveorDead")


Multivariateclinicalcovariates <- TMACLINmerge[,clincovariates]

#Merge the optimal cutpoints with clinical covariates

df_list <- list(COXphoptimalcutoffs, Multivariateclinicalcovariates)

COXphoptimalcutoffsClinMerge <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list)

#simplify values
COXphoptimalcutoffsClinMerge$pT[COXphoptimalcutoffsClinMerge$pT %in% c("1", "1a", "1b") ] <- "0-1"

COXphoptimalcutoffsClinMerge$pT[COXphoptimalcutoffsClinMerge$pT %in% c("2", "2", "3", "4", "4a") ] <- "2-4"

COXphoptimalcutoffsClinMerge$MandardScore[COXphoptimalcutoffsClinMerge$MandardScore %in% c("TRG1", "TRG2") ] <- "TRG1-2"
COXphoptimalcutoffsClinMerge$MandardScore[COXphoptimalcutoffsClinMerge$MandardScore %in% c("TRG3", "TRG4", "TRG5") ] <- "TRG3-5"
COXphoptimalcutoffsClinMerge$MandardScore[COXphoptimalcutoffsClinMerge$MandardScore %in% c("9", "unknown") ] <- "-"
COXphoptimalcutoffsClinMerge$Perineural_Invasion[COXphoptimalcutoffsClinMerge$Perineural_Invasion %in% "2" ] <- "N"
COXphoptimalcutoffsClinMerge$Perineural_Invasion[COXphoptimalcutoffsClinMerge$Perineural_Invasion %in% "1" ] <- "Y"
COXphoptimalcutoffsClinMerge$pN[COXphoptimalcutoffsClinMerge$pN %in% c("1", "2", "3") ] <- "1-3"
COXphoptimalcutoffsClinMerge$pN[COXphoptimalcutoffsClinMerge$pN %in% "0" ] <- "0"
COXphoptimalcutoffsClinMerge$TumourGradingDifferentiation[COXphoptimalcutoffsClinMerge$TumourGradingDifferentiation %in% c("G1well", "G2moderate") ] <- "G1-2"
COXphoptimalcutoffsClinMerge$TumourGradingDifferentiation[COXphoptimalcutoffsClinMerge$TumourGradingDifferentiation %in% c("G3poor", "G4anaplastic") ] <- "G3-4"
COXphoptimalcutoffsClinMerge$TumourGradingDifferentiation[COXphoptimalcutoffsClinMerge$TumourGradingDifferentiation %in% "unknown" ] <- "-"
COXphoptimalcutoffsClinMerge$pM[COXphoptimalcutoffsClinMerge$pM %in% "X" ] <- "-"

COXphoptimalcutoffsClinMerge <- COXphoptimalcutoffsClinMerge %>%
  mutate(across(where(is.character), ~na_if(., "-")))
covariates <- c(    "Age",
                    "Sex",
                    "pT",
                    "pN",
                    "pM")


#Univariate clinical model
univ_formulas <- sapply(covariates, function(x) as.formula(paste('Surv(OS.months, Survival.status)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = COXphoptimalcutoffsClinMerge, na.action = na.omit)})

forest_model(model_list = univ_models,
             covariates = covariates,
             merge_models = T,
             recalculate_height = 15,
             recalculate_width = 3)
#Multivariate clinical model
library(finalfit)
dependent_os  <- "Surv(OS.months, Survival.status)"
COXphoptimalcutoffsClinMerge$CSDE1Hscore
explanatory   <- c(    "Age",
                       "Sex",
                       "pT",
                       "pN",
                       "pM",
                      "HLA_ABCHscore",
                      "HLA_C2Hscore",
                      "HLA_EHscore",
                      "TAP1Hscore",
                      "CSDE1Hscore")
#Relevel if needed
COXphoptimalcutoffsClinMerge$HLA_C2Hscore <- as.factor(COXphoptimalcutoffsClinMerge$HLA_C2Hscore)
COXphoptimalcutoffsClinMerge$HLA_C2Hscore = relevel(COXphoptimalcutoffsClinMerge$HLA_C2Hscore, ref = "low")

COXphoptimalcutoffsClinMerge %>% 
  finalfit(dependent_os, explanatory, metrics = TRUE)

COXphoptimalcutoffsClinMerge %>% 
  finalfit(dependent_os, explanatory, add_dependent_label = FALSE) %>% 
  rename("Overall survival" = label) %>% 
  rename(" " = levels) %>% 
  rename("  " = all)

COXphoptimalcutoffsClinMerge %>% 
  hr_plot(dependent_os, explanatory)

# Heatmap the Stain/count values




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
TMAdemoClin2$Sex
TMAdemoClin2 %>% 
  group_by(Sex) %>% 
  tally()

#median survival

TMAOAC.surv <- as.data.frame(summary((Surv(time = TMACLINmerge$OS.months, event = TMACLINmerge$Survival.status))))
mediansurv <- `TMA datasheet final referenced.csv`
TMAOAC.surv <- as.data.frame(summary((Surv(time = mediansurv$OS.survival.time, event = mediansurv$Survival.status))))


#Disease free survival
setwd("C:/Users/wp1g19/Documents/TMA analysis/R full analysis wd/Maximal concordant survival analysis/DFS")
################# CD3 disease free survival

#make the optimal cut
CD3countpermm2.surv_count.cut.dfs <- surv_cutpoint(
  TMACLINmerge,
  time = "Disease.free.survival.months",
  event = "Reocurrence.status",
  variables = "CD3countpermm2", minprop = 0.45, maxprop = 0.65
)

cd3cutsum.dfs <- as.data.frame(summary(CD3countpermm2.surv_count.cut.dfs))

#Plot the cutpoint 

# palette = "npg" (nature publishing group), see ?ggpubr::ggpar
plot(CD3countpermm2.surv_count.cut.dfs, "CD3countpermm2", palette = "npg", bins = 166)
dev.copy(png,filename="CD3.surv_count.cut.dfs.png");
while (dev.cur()>1) dev.off()

#categorise the cutpoint 
CD3.surv_count.cat.dfs <- surv_categorize(CD3countpermm2.surv_count.cut.dfs) 

library(survival)
library(RTCGA)
fit <- survfit(Surv(Disease.free.survival.months, Reocurrence.status) ~ CD3countpermm2,
               data = CD3.surv_count.cat.dfs)
# extract p value
cd3.dfspvalue <- surv_pvalue(fit)$pval.txt
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,160), 
  xlab = ("Months"),# present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 20,    # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = TRUE # show bars instead of names in text annotations
  # in legend of risk table
)  
#save the image file
dev.copy(png,filename="CD3OSKmDFS.png");
while (dev.cur()>1) dev.off()


################# CD4 disease free survival

#make the optimal cut
CD4countpermm2.surv_count.cut.dfs <- surv_cutpoint(
  TMACLINmerge,
  time = "Disease.free.survival.months",
  event = "Reocurrence.status",
  variables = "CD4countpermm2", minprop = 0.5, maxprop = 0.6
)

CD4cutsum.dfs <- as.data.frame(summary(CD4countpermm2.surv_count.cut.dfs))

#Plot the cutpoint 

# palette = "npg" (nature publishing group), see ?ggpubr::ggpar
plot(CD4countpermm2.surv_count.cut.dfs, "CD4countpermm2", palette = "npg", bins = 166)
dev.copy(png,filename="CD4.surv_count.cut.dfs.png");
while (dev.cur()>1) dev.off()

#categorise the cutpoint 
CD4.surv_count.cat.dfs <- surv_categorize(CD4countpermm2.surv_count.cut.dfs) 

library(survival)
library(RTCGA)
fit <- survfit(Surv(Disease.free.survival.months, Reocurrence.status) ~ CD4countpermm2,
               data = CD4.surv_count.cat.dfs)
# extract p value
CD4.dfspvalue <- surv_pvalue(fit)$pval.txt
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,160), 
  xlab = ("Months"),# present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 20,    # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = TRUE # show bars instead of names in text annotations
  # in legend of risk table
)  
#save the image file
dev.copy(png,filename="CD4OSKmDFS.png");
while (dev.cur()>1) dev.off()

################# CD8 disease free survival

#make the optimal cut
CD8countpermm2.surv_count.cut.dfs <- surv_cutpoint(
  TMACLINmerge,
  time = "Disease.free.survival.months",
  event = "Reocurrence.status",
  variables = "CD8countpermm2", minprop = 0.25, maxprop = 0.5
)

CD8cutsum.dfs <- as.data.frame(summary(CD8countpermm2.surv_count.cut.dfs))

#Plot the cutpoint 

# palette = "npg" (nature publishing group), see ?ggpubr::ggpar
plot(CD8countpermm2.surv_count.cut.dfs, "CD8countpermm2", palette = "npg", bins = 166)
dev.copy(png,filename="CD8.surv_count.cut.dfs.png");
while (dev.cur()>1) dev.off()

#categorise the cutpoint 
CD8.surv_count.cat.dfs <- surv_categorize(CD8countpermm2.surv_count.cut.dfs) 

library(survival)
library(RTCGA)
fit <- survfit(Surv(Disease.free.survival.months, Reocurrence.status) ~ CD8countpermm2,
               data = CD8.surv_count.cat.dfs)
# extract p value
CD8.dfspvalue <- surv_pvalue(fit)$pval.txt
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,160), 
  xlab = ("Months"),# present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 20,    # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = TRUE # show bars instead of names in text annotations
  # in legend of risk table
)  
#save the image file
dev.copy(png,filename="CD8OSKmDFS.png");
while (dev.cur()>1) dev.off()

################# FOXP3 disease free survival

#make the optimal cut
FOXP3countpermm2.surv_count.cut.dfs <- surv_cutpoint(
  TMACLINmerge,
  time = "Disease.free.survival.months",
  event = "Reocurrence.status",
  variables = "FOXP3countpermm2", minprop = 0.1, maxprop = 0.9
)

FOXP3cutsum.dfs <- as.data.frame(summary(FOXP3countpermm2.surv_count.cut.dfs))

#Plot the cutpoint 

# palette = "npg" (nature publishing group), see ?ggpubr::ggpar
plot(FOXP3countpermm2.surv_count.cut.dfs, "FOXP3countpermm2", palette = "npg", bins = 166)
dev.copy(png,filename="FOXP3.surv_count.cut.dfs.png");
while (dev.cur()>1) dev.off()

#categorise the cutpoint 
FOXP3.surv_count.cat.dfs <- surv_categorize(FOXP3countpermm2.surv_count.cut.dfs) 

library(survival)
library(RTCGA)
fit <- survfit(Surv(Disease.free.survival.months, Reocurrence.status) ~ FOXP3countpermm2,
               data = FOXP3.surv_count.cat.dfs)
# extract p value
FOXP3.dfspvalue <- surv_pvalue(fit)$pval.txt
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,160), 
  xlab = ("Months"),# present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 20,    # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = TRUE # show bars instead of names in text annotations
  # in legend of risk table
)  
#save the image file
dev.copy(png,filename="FOXP3OSKmDFS.png");
while (dev.cur()>1) dev.off()

########### HLA-ABC percent disease free survival

#make the optimal cut
HLA_ABCpercentagescore.surv_count.cut.dfs <- surv_cutpoint(
  TMACLINmerge,
  time = "Disease.free.survival.months",
  event = "Reocurrence.status",
  variables = "HLA_ABCpercentagescore", minprop = 0.2, maxprop = 0.6
)

HLA_ABC_perccutsumDFS <- as.data.frame(summary(HLA_ABCpercentagescore.surv_count.cut.dfs))

#Plot the cutpoint 

# palette = "npg" (nature publishing group), see ?ggpubr::ggpar
plot(HLA_ABCpercentagescore.surv_count.cut.dfs, "HLA_ABCpercentagescore", palette = "npg", bins = 166)
dev.copy(png,filename="HLA_ABCpercentagescore.surv_count.cut.dfs.png");
while (dev.cur()>1) dev.off()

#categorise the cutpoint 
HLA_ABCpercentagescore.surv_count.cat.dfs <- surv_categorize(HLA_ABCpercentagescore.surv_count.cut.dfs) 

library(survival)
library(RTCGA)
fit <- survfit(Surv(Disease.free.survival.months, Reocurrence.status) ~ HLA_ABCpercentagescore,
               data = HLA_ABCpercentagescore.surv_count.cat.dfs)
# extract p value
HLAABCperpvalueDFS <- surv_pvalue(fit)$pval.txt
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,180), 
  xlab = ("Months"),# present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 20,    # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = TRUE # show bars instead of names in text annotations
  # in legend of risk table
)  
#save the image file
dev.copy(png,filename="HLA_ABCpercentagescoreOSKmDFS.png");
while (dev.cur()>1) dev.off()


############# HLA_ABC H score disease free survival
#make the optimal cut
HLA_ABCHscore.surv_count.cut.dfs <- surv_cutpoint(
  TMACLINmerge,
  time = "Disease.free.survival.months",
  event = "Reocurrence.status",
  variables = "HLA_ABCHscore", minprop = 0.3, maxprop = 0.7
)

HLA_ABC_HScutsumDFS <- as.data.frame(summary(HLA_ABCHscore.surv_count.cut.dfs))

#Plot the cutpoint 

# palette = "npg" (nature publishing group), see ?ggpubr::ggpar
plot(HLA_ABCHscore.surv_count.cut.dfs, "HLA_ABCHscore", palette = "npg", bins = 166)
dev.copy(png,filename="HLA_ABCHscore.surv_count.cut.dfs.png");
while (dev.cur()>1) dev.off()

#categorise the cutpoint 
HLA_ABCHscore.surv_count.cat.dfs <- surv_categorize(HLA_ABCHscore.surv_count.cut.dfs) 

library(survival)
library(RTCGA)
fit <- survfit(Surv(Disease.free.survival.months, Reocurrence.status) ~ HLA_ABCHscore,
               data = HLA_ABCHscore.surv_count.cat.dfs)
# extract p value
HLAABCHSpvalueDFS <- surv_pvalue(fit)$pval.txt
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,180), 
  xlab = ("Months"),# present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 20,    # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = TRUE # show bars instead of names in text annotations
  # in legend of risk table
)  
#save the image file
dev.copy(png,filename="HLA_ABCHscoreOSKmDFS.png");
while (dev.cur()>1) dev.off()



############# HLA_C2 percentage disease free survival
#make the optimal cut
HLA_C2percentagescore.surv_count.cut.dfs <- surv_cutpoint(
  TMACLINmerge,
  time = "Disease.free.survival.months",
  event = "Reocurrence.status",
  variables = "HLA_C2percentagescore", minprop = 0.25, maxprop = 0.75
)

HLA_C2_perccutsumDFS <- as.data.frame(summary(HLA_C2percentagescore.surv_count.cut.dfs))

#Plot the cutpoint 

# palette = "npg" (nature publishing group), see ?ggpubr::ggpar
plot(HLA_C2percentagescore.surv_count.cut.dfs, "HLA_C2percentagescore", palette = "npg", bins = 166)
dev.copy(png,filename="HLA_C2percentagescore.surv_count.cut.dfs.png");
while (dev.cur()>1) dev.off()

#categorise the cutpoint 
HLA_C2percentagescore.surv_count.cat.dfs <- surv_categorize(HLA_C2percentagescore.surv_count.cut.dfs) 

library(survival)
library(RTCGA)
fit <- survfit(Surv(Disease.free.survival.months, Reocurrence.status) ~ HLA_C2percentagescore,
               data = HLA_C2percentagescore.surv_count.cat.dfs)
# extract p value
HLAC2perpvalueDFS <- surv_pvalue(fit)$pval.txt
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,180), 
  xlab = ("Months"),# present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 20,    # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = TRUE # show bars instead of names in text annotations
  # in legend of risk table
)  
#save the image file
dev.copy(png,filename="HLA_C2percentagescoreOSKmDFS.png");
while (dev.cur()>1) dev.off()

############# HLA_C2 H score disease free survival
#make the optimal cut
HLA_C2Hscore.surv_count.cut.dfs <- surv_cutpoint(
  TMACLINmerge,
  time = "Disease.free.survival.months",
  event = "Reocurrence.status",
  variables = "HLA_C2Hscore", minprop = 0.1, maxprop = 0.9
)

HLA_C2_HScutsumDFS <- as.data.frame(summary(HLA_C2Hscore.surv_count.cut.dfs))

#Plot the cutpoint 

# palette = "npg" (nature publishing group), see ?ggpubr::ggpar
plot(HLA_C2Hscore.surv_count.cut.dfs, "HLA_C2Hscore", palette = "npg", bins = 166)
dev.copy(png,filename="HLA_C2Hscore.surv_count.cut.dfs.png");
while (dev.cur()>1) dev.off()

#categorise the cutpoint 
HLA_C2Hscore.surv_count.cat.dfs <- surv_categorize(HLA_C2Hscore.surv_count.cut.dfs) 

library(survival)
library(RTCGA)
fit <- survfit(Surv(Disease.free.survival.months, Reocurrence.status) ~ HLA_C2Hscore,
               data = HLA_C2Hscore.surv_count.cat.dfs)
# extract p value
HLAC2HSpvalueDFS <- surv_pvalue(fit)$pval.txt
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,180), 
  xlab = ("Months"),# present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 20,    # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = TRUE # show bars instead of names in text annotations
  # in legend of risk table
)  
#save the image file
dev.copy(png,filename="HLA_C2HscoreOSKmDFS.png");
while (dev.cur()>1) dev.off()

############# HLA_Epercentagescore disease free survival
#make the optimal cut
HLA_Epercentagescore.surv_count.cut.dfs <- surv_cutpoint(
  TMACLINmerge,
  time = "Disease.free.survival.months",
  event = "Reocurrence.status",
  variables = "HLA_Epercentagescore", minprop = 0.5, maxprop = 0.8
)

HLA_E_perccutsumDFS <- as.data.frame(summary(HLA_Epercentagescore.surv_count.cut.dfs))

#Plot the cutpoint 

# palette = "npg" (nature publishing group), see ?ggpubr::ggpar
plot(HLA_Epercentagescore.surv_count.cut.dfs, "HLA_Epercentagescore", palette = "npg", bins = 166)
dev.copy(png,filename="HLA_Epercentagescore.surv_count.cut.dfs.png");
while (dev.cur()>1) dev.off()

#categorise the cutpoint 
HLA_Epercentagescore.surv_count.cat.dfs <- surv_categorize(HLA_Epercentagescore.surv_count.cut.dfs) 

library(survival)
library(RTCGA)
fit <- survfit(Surv(Disease.free.survival.months, Reocurrence.status) ~ HLA_Epercentagescore,
               data = HLA_Epercentagescore.surv_count.cat.dfs)
# extract p value
HLAEperpvalueDFS <- surv_pvalue(fit)$pval.txt
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,180), 
  xlab = ("Months"),# present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 20,    # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = TRUE # show bars instead of names in text annotations
  # in legend of risk table
)  
#save the image file
dev.copy(png,filename="HLA_EpercentagescoreOSKmDFS.png");
while (dev.cur()>1) dev.off()

############# HLA_EHscore survival
#make the optimal cut
HLA_EHscore.surv_count.cut.dfs <- surv_cutpoint(
  TMACLINmerge,
  time = "Disease.free.survival.months",
  event = "Reocurrence.status",
  variables = "HLA_EHscore", minprop = 0.2, maxprop = 0.75
)

HLA_E_HScutsumDFS <- as.data.frame(summary(HLA_EHscore.surv_count.cut.dfs))

#Plot the cutpoint 

# palette = "npg" (nature publishing group), see ?ggpubr::ggpar
plot(HLA_EHscore.surv_count.cut.dfs, "HLA_EHscore", palette = "npg", bins = 166)
dev.copy(png,filename="HLA_EHscore.surv_count.cut.dfs.png");
while (dev.cur()>1) dev.off()

#categorise the cutpoint 
HLA_EHscore.surv_count.cat.dfs <- surv_categorize(HLA_EHscore.surv_count.cut.dfs) 

library(survival)
library(RTCGA)
fit <- survfit(Surv(Disease.free.survival.months, Reocurrence.status) ~ HLA_EHscore,
               data = HLA_EHscore.surv_count.cat.dfs)
# extract p value
HLAEHSpvalueDFS <- surv_pvalue(fit)$pval.txt
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,180), 
  xlab = ("Months"),# present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 20,    # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = TRUE # show bars instead of names in text annotations
  # in legend of risk table
)  
#save the image file
dev.copy(png,filename="HLA_EHscoreOSKmDFS.png");
while (dev.cur()>1) dev.off()

############# TAP1percentagescore disease free survival
#make the optimal cut
TAP1percentagescore.surv_count.cut.dfs <- surv_cutpoint(
  TMACLINmerge,
  time = "Disease.free.survival.months",
  event = "Reocurrence.status",
  variables = "TAP1percentagescore", minprop = 0.3, maxprop = 0.7
)

TAP1_perccutsumDFS <- as.data.frame(summary(TAP1percentagescore.surv_count.cut.dfs))

#Plot the cutpoint 

# palette = "npg" (nature publishing group), see ?ggpubr::ggpar
plot(TAP1percentagescore.surv_count.cut.dfs, "TAP1percentagescore", palette = "npg", bins = 166)
dev.copy(png,filename="TAP1percentagescore.surv_count.cut.dfs.png");
while (dev.cur()>1) dev.off()

#categorise the cutpoint 
TAP1percentagescore.surv_count.cat.dfs <- surv_categorize(TAP1percentagescore.surv_count.cut.dfs) 

library(survival)
library(RTCGA)
fit <- survfit(Surv(Disease.free.survival.months, Reocurrence.status) ~ TAP1percentagescore,
               data = TAP1percentagescore.surv_count.cat.dfs)
# extract p value
TAP1perpvalueDFS <- surv_pvalue(fit)$pval.txt
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,180), 
  xlab = ("Months"),# present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 20,    # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = TRUE # show bars instead of names in text annotations
  # in legend of risk table
)  
#save the image file
dev.copy(png,filename="TAP1percentagescoreOSKmDFS.png");
while (dev.cur()>1) dev.off()

############# TAP1Hscore disease free survival
#make the optimal cut
TAP1Hscore.surv_count.cut.dfs <- surv_cutpoint(
  TMACLINmerge,
  time = "Disease.free.survival.months",
  event = "Reocurrence.status",
  variables = "TAP1Hscore", minprop = 0.25, maxprop = 0.7
)

TAP1_HScutsumDFS <- as.data.frame(summary(TAP1Hscore.surv_count.cut.dfs))

#Plot the cutpoint 

# palette = "npg" (nature publishing group), see ?ggpubr::ggpar
plot(TAP1Hscore.surv_count.cut.dfs, "TAP1Hscore", palette = "npg", bins = 166)
dev.copy(png,filename="TAP1Hscore.surv_count.cut.dfs.png");
while (dev.cur()>1) dev.off()

#categorise the cutpoint 
TAP1Hscore.surv_count.cat.dfs <- surv_categorize(TAP1Hscore.surv_count.cut.dfs) 

library(survival)
library(RTCGA)
fit <- survfit(Surv(Disease.free.survival.months, Reocurrence.status) ~ TAP1Hscore,
               data = TAP1Hscore.surv_count.cat.dfs)
# extract p value
TAP1HSpvalueDFS <- surv_pvalue(fit)$pval.txt

ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,180), 
  xlab = ("Months"),# present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 20,    # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = TRUE # show bars instead of names in text annotations
  # in legend of risk table
)  
#save the image file
dev.copy(png,filename="TAP1HscoreOSKmDFS.png");
while (dev.cur()>1) dev.off()

cutpointsums.dfs <- rbind(cd3cutsum.dfs,
                      CD4cutsum.dfs,
                      CD8cutsum.dfs,
                      FOXP3cutsum.dfs,
                      HLA_ABC_perccutsumDFS,
                      HLA_ABC_HScutsumDFS,
                      HLA_C2_perccutsumDFS,
                      HLA_C2_HScutsumDFS,
                      HLA_E_perccutsumDFS,
                      HLA_E_HScutsumDFS,
                      TAP1_perccutsumDFS,
                      TAP1_HScutsumDFS)
pvalue.dfs <- as.data.frame(rbind(cd3.dfspvalue,
                               CD4.dfspvalue,
                               CD8.dfspvalue,
                               FOXP3.dfspvalue,
                               HLAABCperpvalueDFS,
                               HLAABCHSpvalueDFS,
                               HLAC2perpvalueDFS,
                               HLAC2HSpvalueDFS,
                               HLAEperpvalueDFS,
                               HLAEHSpvalueDFS,
                               TAP1perpvalueDFS,
                               TAP1HSpvalueDFS))
cutpointpvaluesums.dfs <- as.data.frame(cutpointsums.dfs)
cutpointpvaluesums.dfs$pvalueDFS <- pvalue.dfs
cutpointpvaluesums.OS.dfs <- merge(cutpointpvaluesums, cutpointpvaluesums.dfs, by = 0)


#Cancer specific survival

library(survival)

survival_object = TMACLINmerge %$% 
  Surv(OS.months, CSS.status)

TMACLINmerge$CSS.status
my_survfit = survfit(survival_object ~ 1, data = TMACLINmerge)
dependent_os = "Surv(OS.months)"
explanatory = c("CSS.status")

TMACLINmerge %>% 
  surv_plot(dependent_os,
            explanatory,
            pval = TRUE,
            conf.int = TRUE,
            break.time.by = 20,  xlab = ("Months"), 
            ggtheme = theme_RTCGA(),
            surv.median.line = "hv",
            linetype = 1)
#save the image file
dev.copy(png,filename="CSS.png");
while (dev.cur()>1) dev.off()

TMACLINmerge$CSS.status.12 <- TMACLINmerge$CSS.status
# TMACLINmerge$CSS.status.12 <- TMACLINmerge$CSS.status <- gsub('3','0',TMACLINmerge$CSS.status)
TMACLINmerge <- TMACLINmerge[TMACLINmerge$CSS.status.12 != 3, ]
TMACLINmerge$CSS.status.12 <- as.integer(TMACLINmerge$CSS.status.12)
setwd("C:/Users/wp1g19/Documents/TMA analysis/R full analysis wd/Maximal concordant survival analysis/CSS")
write.csv(TMACLINmerge, file = "C:/Users/wp1g19/Documents/TMA analysis/R full analysis wd/images for thesis/Quantile_cut/CSS remove death by other/TMAclinmergeKMplotterDFS.csv")
################# CD3 Cancer specific survival

#make the optimal cut
CD3countpermm2.surv_count.cut.css <- surv_cutpoint(
  TMACLINmerge,
  time = "OS.months",
  event = "CSS.status.12",
  variables = "CD3countpermm2", minprop = 0.45, maxprop = 0.65
)

cd3cutsum.CSS <- as.data.frame(summary(CD3countpermm2.surv_count.cut.css))

#Plot the cutpoint 

# palette = "npg" (nature publishing group), see ?ggpubr::ggpar
plot(CD3countpermm2.surv_count.cut.css, "CD3countpermm2", palette = "npg", bins = 166)
dev.copy(png,filename="CD3.surv_count.cut.css.png");
while (dev.cur()>1) dev.off()

#categorise the cutpoint 
CD3.surv_count.cat.css <- surv_categorize(CD3countpermm2.surv_count.cut.css) 

library(survival)
library(RTCGA)

fit <- survfit(Surv(time = OS.months,event = CSS.status.12) ~ CD3countpermm2, data = CD3.surv_count.cat.css)

# extract p value
cd3.csspvalue <- surv_pvalue(fit)$pval.txt
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,160), 
  xlab = ("Months"),# present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 20,    # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = TRUE # show bars instead of names in text annotations
  # in legend of risk table
)  
#save the image file
dev.copy(png,filename="CD3OSKmCSS.png");
while (dev.cur()>1) dev.off()


################# CD4 Cancer specific survival

#make the optimal cut
CD4countpermm2.surv_count.cut.css <- surv_cutpoint(
  TMACLINmerge,
  time = "OS.months",
  event = "CSS.status.12",
  variables = "CD4countpermm2", minprop = 0.3, maxprop = 0.7
)

CD4cutsum.CSS <- as.data.frame(summary(CD4countpermm2.surv_count.cut.css))

#Plot the cutpoint 

# palette = "npg" (nature publishing group), see ?ggpubr::ggpar
plot(CD4countpermm2.surv_count.cut.css, "CD4countpermm2", palette = "npg", bins = 166)
dev.copy(png,filename="CD4.surv_count.cut.css.png");
while (dev.cur()>1) dev.off()

#categorise the cutpoint 
CD4.surv_count.cat.css <- surv_categorize(CD4countpermm2.surv_count.cut.css) 

library(survival)
library(RTCGA)

fit <- survfit(Surv(time = OS.months,event = CSS.status.12) ~ CD4countpermm2, data = CD4.surv_count.cat.css)

# extract p value
CD4.csspvalue <- surv_pvalue(fit)$pval.txt
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,160), 
  xlab = ("Months"),# present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 20,    # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = TRUE # show bars instead of names in text annotations
  # in legend of risk table
)  
#save the image file
dev.copy(png,filename="CD4OSKmCSS.png");
while (dev.cur()>1) dev.off()

################# CD8 Cancer specific survival

#make the optimal cut
CD8countpermm2.surv_count.cut.css <- surv_cutpoint(
  TMACLINmerge,
  time = "OS.months",
  event = "CSS.status.12",
  variables = "CD8countpermm2", minprop = 0.2, maxprop = 0.8
)

CD8cutsum.CSS <- as.data.frame(summary(CD8countpermm2.surv_count.cut.css))

#Plot the cutpoint 

# palette = "npg" (nature publishing group), see ?ggpubr::ggpar
plot(CD8countpermm2.surv_count.cut.css, "CD8countpermm2", palette = "npg", bins = 166)
dev.copy(png,filename="CD8.surv_count.cut.css.png");
while (dev.cur()>1) dev.off()

#categorise the cutpoint 
CD8.surv_count.cat.css <- surv_categorize(CD8countpermm2.surv_count.cut.css) 

library(survival)
library(RTCGA)

fit <- survfit(Surv(time = OS.months,event = CSS.status.12) ~ CD8countpermm2, data = CD8.surv_count.cat.css)

# extract p value
CD8.csspvalue <- surv_pvalue(fit)$pval.txt
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,160), 
  xlab = ("Months"),# present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 20,    # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = TRUE # show bars instead of names in text annotations
  # in legend of risk table
)  
#save the image file
dev.copy(png,filename="CD8OSKmCSS.png");
while (dev.cur()>1) dev.off()

################# FOXP3 Cancer specific survival

#make the optimal cut
FOXP3countpermm2.surv_count.cut.css <- surv_cutpoint(
  TMACLINmerge,
  time = "OS.months",
  event = "CSS.status.12",
  variables = "FOXP3countpermm2", minprop = 0.1, maxprop = 0.9
)

FOXP3cutsum.CSS <- as.data.frame(summary(FOXP3countpermm2.surv_count.cut.css))

#Plot the cutpoint 

# palette = "npg" (nature publishing group), see ?ggpubr::ggpar
plot(FOXP3countpermm2.surv_count.cut.css, "FOXP3countpermm2", palette = "npg", bins = 166)
dev.copy(png,filename="FOXP3.surv_count.cut.css.png");
while (dev.cur()>1) dev.off()

#categorise the cutpoint 
FOXP3.surv_count.cat.css <- surv_categorize(FOXP3countpermm2.surv_count.cut.css) 

library(survival)
library(RTCGA)

fit <- survfit(Surv(time = OS.months,event = CSS.status.12) ~ FOXP3countpermm2, data = FOXP3.surv_count.cat.css)

# extract p value
FOXP3.csspvalue <- surv_pvalue(fit)$pval.txt
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,160), 
  xlab = ("Months"),# present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 20,    # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = TRUE # show bars instead of names in text annotations
  # in legend of risk table
)  
#save the image file
dev.copy(png,filename="FOXP3OSKmCSS.png");
while (dev.cur()>1) dev.off()



########### HLA-ABC percent cancer specific survival

#make the optimal cut
HLA_ABCpercentagescore.surv_count.cut.css <- surv_cutpoint(
  TMACLINmerge,
  time = "OS.months",
  event = "CSS.status.12",
  variables = "HLA_ABCpercentagescore", minprop = 0.25, maxprop = 0.6
)

HLA_ABC_perccutsumCSS <- as.data.frame(summary(HLA_ABCpercentagescore.surv_count.cut.css))

#Plot the cutpoint 

# palette = "npg" (nature publishing group), see ?ggpubr::ggpar
plot(HLA_ABCpercentagescore.surv_count.cut.css, "HLA_ABCpercentagescore", palette = "npg", bins = 166)
dev.copy(png,filename="HLA_ABCpercentagescore.surv_count.cut.css.png");
while (dev.cur()>1) dev.off()

#categorise the cutpoint 
HLA_ABCpercentagescore.surv_count.cat.css <- surv_categorize(HLA_ABCpercentagescore.surv_count.cut.css) 

library(survival)
library(RTCGA)
fit <- survfit(Surv(OS.months, CSS.status.12) ~ HLA_ABCpercentagescore,
               data = HLA_ABCpercentagescore.surv_count.cat.css)
# extract p value
HLAABCperpvalueCSS <- surv_pvalue(fit)$pval.txt
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,180), 
  xlab = ("Months"),# present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 20,    # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = TRUE # show bars instead of names in text annotations
  # in legend of risk table
)  
#save the image file
dev.copy(png,filename="HLA_ABCpercentagescoreCSSKM.png");
while (dev.cur()>1) dev.off()


############# HLA_ABC H score cancer specific survival
#make the optimal cut
HLA_ABCHscore.surv_count.cut.css <- surv_cutpoint(
  TMACLINmerge,
  time = "OS.months",
  event = "CSS.status.12",
  variables = "HLA_ABCHscore", minprop = 0.4, maxprop = 0.8
)

HLA_ABC_HScutsumCSS <- as.data.frame(summary(HLA_ABCHscore.surv_count.cut.css))

#Plot the cutpoint 

# palette = "npg" (nature publishing group), see ?ggpubr::ggpar
plot(HLA_ABCHscore.surv_count.cut.css, "HLA_ABCHscore", palette = "npg", bins = 166)
dev.copy(png,filename="HLA_ABCHscore.surv_count.cut.css.png");
while (dev.cur()>1) dev.off()

#categorise the cutpoint 
HLA_ABCHscore.surv_count.cat.css <- surv_categorize(HLA_ABCHscore.surv_count.cut.css) 

library(survival)
library(RTCGA)
fit <- survfit(Surv(OS.months, CSS.status.12) ~ HLA_ABCHscore,
               data = HLA_ABCHscore.surv_count.cat.css)
# extract p value
HLAABCHSpvalueCSS <- surv_pvalue(fit)$pval.txt
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,60), 
  xlab = ("Months"),# present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 20,    # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), font.x = 25, font.y = 25, font.legend = 25, font.caption = 25,  font.tickslab = 12, # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = TRUE, surv.median.line = "hv" # show bars instead of names in text annotations
  # in legend of risk table
)  + 
  guides(colour = guide_legend(nrow = 2))
#save the image file
dev.copy(png,filename="HLA_ABCHscoreCSSKM.png");
while (dev.cur()>1) dev.off()



############# HLA_C2 percent cancer specific survival
#make the optimal cut
HLA_C2percentagescore.surv_count.cut.css <- surv_cutpoint(
  TMACLINmerge,
  time = "OS.months",
  event = "CSS.status.12",
  variables = "HLA_C2percentagescore", minprop = 0.32, maxprop = 0.7
)

HLA_C2_perccutsumCSS <- as.data.frame(summary(HLA_C2percentagescore.surv_count.cut.css))

#Plot the cutpoint 

# palette = "npg" (nature publishing group), see ?ggpubr::ggpar
plot(HLA_C2percentagescore.surv_count.cut.css, "HLA_C2percentagescore", palette = "npg", bins = 166)
dev.copy(png,filename="HLA_C2percentagescore.surv_count.cut.css.png");
while (dev.cur()>1) dev.off()

#categorise the cutpoint 
HLA_C2percentagescore.surv_count.cat.css <- surv_categorize(HLA_C2percentagescore.surv_count.cut.css) 

library(survival)
library(RTCGA)
fit <- survfit(Surv(OS.months, CSS.status.12) ~ HLA_C2percentagescore,
               data = HLA_C2percentagescore.surv_count.cat.css)
# extract p value
HLAC2perpvalueCSS <- surv_pvalue(fit)$pval.txt
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,180), 
  xlab = ("Months"),# present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 20,    # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = TRUE # show bars instead of names in text annotations
  # in legend of risk table
)  
#save the image file
dev.copy(png,filename="HLA_C2percentagescoreCSSKM.png");
while (dev.cur()>1) dev.off()

############# HLA_C2 H score cancer specific survival
#make the optimal cut
HLA_C2Hscore.surv_count.cut.css <- surv_cutpoint(
  TMACLINmerge,
  time = "OS.months",
  event = "CSS.status.12",
  variables = "HLA_C2Hscore", minprop = 0.1, maxprop = 0.9
)

HLA_C2_HScutsumCSS <- as.data.frame(summary(HLA_C2Hscore.surv_count.cut.css))

#Plot the cutpoint 

# palette = "npg" (nature publishing group), see ?ggpubr::ggpar
plot(HLA_C2Hscore.surv_count.cut.css, "HLA_C2Hscore", palette = "npg", bins = 166)
dev.copy(png,filename="HLA_C2Hscore.surv_count.cut.css.png");
while (dev.cur()>1) dev.off()

#categorise the cutpoint 
HLA_C2Hscore.surv_count.cat.css <- surv_categorize(HLA_C2Hscore.surv_count.cut.css) 

library(survival)
library(RTCGA)
fit <- survfit(Surv(OS.months, CSS.status.12) ~ HLA_C2Hscore,
               data = HLA_C2Hscore.surv_count.cat.css)
# extract p value
HLAC2HSpvalueCSS <- surv_pvalue(fit)$pval.txt
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,180), 
  xlab = ("Months"),# present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 20,    # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = TRUE # show bars instead of names in text annotations
  # in legend of risk table
)  
#save the image file
dev.copy(png,filename="HLA_C2HscoreCSSKM.png");
while (dev.cur()>1) dev.off()

############# HLA_Epercent cancer specific survival
#make the optimal cut
HLA_Epercentagescore.surv_count.cut.css <- surv_cutpoint(
  TMACLINmerge,
  time = "OS.months",
  event = "CSS.status.12",
  variables = "HLA_Epercentagescore", minprop = 0.4
)

HLA_E_perccutsumCSS <- as.data.frame(summary(HLA_Epercentagescore.surv_count.cut.css))

#Plot the cutpoint 

# palette = "npg" (nature publishing group), see ?ggpubr::ggpar
plot(HLA_Epercentagescore.surv_count.cut.css, "HLA_Epercentagescore", palette = "npg", bins = 166)
dev.copy(png,filename="HLA_Epercentagescore.surv_count.cut.css.png");
while (dev.cur()>1) dev.off()

#categorise the cutpoint 
HLA_Epercentagescore.surv_count.cat.css <- surv_categorize(HLA_Epercentagescore.surv_count.cut.css) 

library(survival)
library(RTCGA)
fit <- survfit(Surv(OS.months, CSS.status.12) ~ HLA_Epercentagescore,
               data = HLA_Epercentagescore.surv_count.cat.css)
# extract p value
HLAEperpvalueCSS <- surv_pvalue(fit)$pval.txt
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,180), 
  xlab = ("Months"),# present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 20,    # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = TRUE # show bars instead of names in text annotations
  # in legend of risk table
)  
#save the image file
dev.copy(png,filename="HLA_EpercentagescoreCSSKM.png");
while (dev.cur()>1) dev.off()

############# HLA_EH score cancer specific survival
#make the optimal cut
HLA_EHscore.surv_count.cut.css <- surv_cutpoint(
  TMACLINmerge,
  time = "OS.months",
  event = "CSS.status.12",
  variables = "HLA_EHscore", minprop = 0.3, maxprop = 0.7
)

HLA_E_HScutsumCSS <- as.data.frame(summary(HLA_EHscore.surv_count.cut.css))

#Plot the cutpoint 

# palette = "npg" (nature publishing group), see ?ggpubr::ggpar
plot(HLA_EHscore.surv_count.cut.css, "HLA_EHscore", palette = "npg", bins = 166)
dev.copy(png,filename="HLA_EHscore.surv_count.cut.css.png");
while (dev.cur()>1) dev.off()

#categorise the cutpoint 
HLA_EHscore.surv_count.cat.css <- surv_categorize(HLA_EHscore.surv_count.cut.css) 

library(survival)
library(RTCGA)
fit <- survfit(Surv(OS.months, CSS.status.12) ~ HLA_EHscore,
               data = HLA_EHscore.surv_count.cat.css)
# extract p value
HLAEHSpvalueCSS <- surv_pvalue(fit)$pval.txt
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,180), 
  xlab = ("Months"),# present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 20,    # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = TRUE # show bars instead of names in text annotations
  # in legend of risk table
)  
#save the image file
dev.copy(png,filename="HLA_EHscoreCSSKM.png");
while (dev.cur()>1) dev.off()

############# TAP1percent cancer specific survival
#make the optimal cut
TAP1percentagescore.surv_count.cut.css <- surv_cutpoint(
  TMACLINmerge,
  time = "OS.months",
  event = "CSS.status.12",
  variables = "TAP1percentagescore", minprop = 0.3, maxprop = 0.7
)

TAP1_perccutsumCSS <- as.data.frame(summary(TAP1percentagescore.surv_count.cut.css))

#Plot the cutpoint 

# palette = "npg" (nature publishing group), see ?ggpubr::ggpar
plot(TAP1percentagescore.surv_count.cut.css, "TAP1percentagescore", palette = "npg", bins = 166)
dev.copy(png,filename="TAP1percentagescore.surv_count.cut.css.png");
while (dev.cur()>1) dev.off()

#categorise the cutpoint 
TAP1percentagescore.surv_count.cat.css <- surv_categorize(TAP1percentagescore.surv_count.cut.css) 

library(survival)
library(RTCGA)
fit <- survfit(Surv(OS.months, CSS.status.12) ~ TAP1percentagescore,
               data = TAP1percentagescore.surv_count.cat.css)
# extract p value
TAP1perpvalueCSS <- surv_pvalue(fit)$pval.txt
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,180), 
  xlab = ("Months"),# present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 20,    # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = TRUE # show bars instead of names in text annotations
  # in legend of risk table
)  
#save the image file
dev.copy(png,filename="TAP1percentagescoreCSSKM.png");
while (dev.cur()>1) dev.off()

############# TAP1Hscore disease free survival
#make the optimal cut
TAP1Hscore.surv_count.cut.css <- surv_cutpoint(
  TMACLINmerge,
  time = "OS.months",
  event = "CSS.status.12",
  variables = "TAP1Hscore", minprop = 0.2, maxprop = 0.6
)

TAP1_HScutsumCSS <- as.data.frame(summary(TAP1Hscore.surv_count.cut.css))

#Plot the cutpoint 

# palette = "npg" (nature publishing group), see ?ggpubr::ggpar
plot(TAP1Hscore.surv_count.cut.css, "TAP1Hscore", palette = "npg", bins = 166)
dev.copy(png,filename="TAP1Hscore.surv_count.cut.css.png");
while (dev.cur()>1) dev.off()

#categorise the cutpoint 
TAP1Hscore.surv_count.cat.css <- surv_categorize(TAP1Hscore.surv_count.cut.css) 

library(survival)
library(RTCGA)
fit <- survfit(Surv(OS.months, CSS.status.12) ~ TAP1Hscore,
               data = TAP1Hscore.surv_count.cat.css)
# extract p value
TAP1HSpvalueCSS <- surv_pvalue(fit)$pval.txt

ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,180), 
  xlab = ("Months"),# present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 20,    # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = TRUE # show bars instead of names in text annotations
  # in legend of risk table
)  
#save the image file
dev.copy(png,filename="TAP1HscoreCSSKM.png");
while (dev.cur()>1) dev.off()



cutpointsums.css <- rbind(cd3cutsum.CSS,
                          CD4cutsum.CSS,
                          CD8cutsum.CSS,
                          FOXP3cutsum.CSS,
                          HLA_ABC_perccutsumCSS,
                          HLA_ABC_HScutsumCSS,
                          HLA_C2_perccutsumCSS,
                          HLA_C2_HScutsumCSS,
                          HLA_E_perccutsumCSS,
                          HLA_E_HScutsumCSS,
                          TAP1_perccutsumCSS,
                          TAP1_HScutsumCSS)
pvalue.css <- as.data.frame(rbind(cd3.csspvalue,
                                  CD4.csspvalue,
                                  CD8.csspvalue,
                                  FOXP3.csspvalue,
                                  HLAABCperpvalueCSS,
                                  HLAABCHSpvalueCSS,
                                  HLAC2perpvalueCSS,
                                  HLAC2HSpvalueCSS,
                                  HLAEperpvalueCSS,
                                  HLAEHSpvalueCSS,
                                  TAP1perpvalueCSS,
                                  TAP1HSpvalueCSS))
cutpointpvaluesums.css <- as.data.frame(cutpointsums.css)
cutpointpvaluesums.css$pvalueCSS <- pvalue.css
cutpointpvaluesums.OS.dfs <- cutpointpvaluesums.OS.dfs[-1]
rownames(cutpointpvaluesums.OS.dfs) <- rownames(cutpointpvaluesums.css)
cutpointpvaluesums.OS.dfs.css <- merge(cutpointpvaluesums.OS.dfs, cutpointpvaluesums.css, by = 0)
#change col/rownames
rnames <- c("CD3permm2",
            "CD4permm2",
            "CD8permm2",
            "FOXP3permm2",
            "HLA_ABC_perc",
            "HLA_ABC_HS",
            "HLA_C2_perc",
            "HLA_C2_HS",
            "HLA_E_perc",
            "HLA_E_HS",
            "TAP1_perc",
            "TAP1_HS")
cnames <- c("cutpointOS",
            "statisticOS",
            "PvalueOS",
            "cutpointDFS",
            "statisticDFS",
            "PvalueDFS",
            "cutpointCSS",
            "statisticCSS",
            "PvalueCSS")
cutpointpvaluesums.OS.dfs.css <- cutpointpvaluesums.OS.dfs.css[-1]
rownames(cutpointpvaluesums.OS.dfs.css) <- rnames
colnames(cutpointpvaluesums.OS.dfs.css) <-cnames
cutpointpvaluesums.OS.dfs.css <- as.data.frame(cutpointpvaluesums.OS.dfs.css)
write.csv(cutpointpvaluesums.OS.dfs.css, file = "C:/Users/wp1g19/Documents/TMA analysis/R full analysis wd/images for thesis/Quantile_cut/cutpointpvaluesums.OS.dfs.css.csv")


#Forest plot coxPH analysis DFS  

clincovariates <- c("Histopathology.Reference.Number",
                    "Age",
                    "Sex",
                    "pT",
                    "pN",
                    "pM",
                    "LymphaticInvasion",
                    "VascularInvasion",
                    "Perineural_Invasion",
                    "TumourGradingDifferentiation",
                    "R0.American.",
                    "MandardScore",
                    "AliveorDead",
                    "Disease.free.survival.months",
                    "Reocurrence.status",
                    "CSS.status.12")

Multivariateclinicalcovariates <- TMACLINmerge[,clincovariates]

#Merge the optimal cutpoints with clinical covariates

df_list <- list(COXphoptimalcutoffs, Multivariateclinicalcovariates)

COXphoptimalcutoffsClinMergeDFSCSS <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list)
write.csv(COXphoptimalcutoffsClinMergeDFSCSS, file = "C:/Users/wp1g19/Documents/TMA analysis/R full analysis wd/images for thesis/Quantile_cut/COXphoptimalcutoffsClinMergeDFSCSS.csv")

setwd()
library(forestmodel)
univ_formulas <- sapply(covariates, function(x) as.formula(paste('Surv(Disease.free.survival.months, Reocurrence.status)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = COXphoptimalcutoffsClinMergeDFSCSS)})

forest_model(model_list = univ_models,
             covariates = covariates,
             merge_models = T,
             recalculate_height = 15,
             return_data = T)
library(grid)
grid.text("Univariate DFS", .3, 0.97, gp=gpar(cex=3, fontsize = 6), check.overlap = T)
dev.copy(png,filename="UnivariateCoxPhForestDFS.png");
while (dev.cur()>1) dev.off()
# Cancer specific univariate coxph

univ_formulas <- sapply(covariates, function(x) as.formula(paste('Surv(OS.months, CSS.status.12)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = COXphoptimalcutoffsClinMergeDFSCSS)})

forest_model(model_list = univ_models,
             covariates = covariates,
             merge_models = T,
             recalculate_height = 15,
             return_data = T)
library(grid)
grid.text("Univariate CSS", .3, 0.97, gp=gpar(cex=3, fontsize = 6), check.overlap = T)
dev.copy(png,filename="UnivariateCoxPhForestCSS.png");
while (dev.cur()>1) dev.off()


#Immune partner APM OS survival
setwd("C:/Users/wp1g19/Documents/TMA analysis/R full analysis wd/Maximal concordant survival analysis/OS_APM_Immunepartners")
#make the optimal cut
optimalcutoffsALL <- COXphoptimalcutoffsClinMerge

COXphoptimalcutoffs$HLA_ABCHscore

#HLA-ABC CD3 OS survival

fit <- survfit(Surv(OS.months, Survival.status) ~ HLA_ABCHscore + CD3countpermm2,
               data = COXphoptimalcutoffs)
# extract p value
HLA_ABCcd3pvalue <- surv_pvalue(fit)$pval.txt
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,180), 
  xlab = ("Months"),# present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 20,    # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = F, # show bars instead of names in text annotations
  # in legend of risk table
  risk.table.height = 0.25,
  legend = "top"
)  + 
  guides(colour = guide_legend(nrow = 2))


#save the image file
dev.copy(png,filename="HLA_ABCCD3OSKm.png");
while (dev.cur()>1) dev.off()



#HLA-ABC CD8 OS survival

fit <- survfit(Surv(OS.months, Survival.status) ~ HLA_ABCHscore + CD8countpermm2,
               data = COXphoptimalcutoffs)
# extract p value
HLA_ABCcd8pvalue <- surv_pvalue(fit)$pval.txt
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,180), 
  xlab = ("Months"),# present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 20,    # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = F, # show bars instead of names in text annotations
  # in legend of risk table
  risk.table.height = 0.3,
  legend = "top"
  )  + 
  guides(colour = guide_legend(nrow = 2))


#save the image file
dev.copy(png,filename="HLA_ABCCD8OSKm.png");
while (dev.cur()>1) dev.off()
#This result of HLA-ABC and CD8 is of high interest Write to excel and remove "HLA_ABClow-CD8high" and "HLA_ABClow-CD8low"
write.csv(COXphoptimalcutoffsClinMerge, file = "C:/Users/wp1g19/Documents/TMA analysis/R full analysis wd/Maximal concordant survival analysis/OS_APM_Immunepartners/OptimalcutpointsALL.csv")
HLA_ABCHighonlyincluded <- read.csv(file = "C:/Users/wp1g19/Documents/TMA analysis/R full analysis wd/Maximal concordant survival analysis/OS_APM_Immunepartners/CutpointsHLA_ABCHighOnly.csv", header = T, row.names = 1)
COXphoptimalcutoffsClinMerge <- read.csv(file = "C:/Users/wp1g19/Documents/TMA analysis/R full analysis wd/Maximal concordant survival analysis/OS_APM_Immunepartners/OptimalcutpointsALL.csv", header = T, row.names = 1)
#HLA-ABC CD8 OS survival HLA high only

HLA_ABCHighonlycuts <- HLA_ABCHighonlyincluded[HLA_ABCHighonlyincluded$HLA_ABCHscore == "high",]
fit <- survfit(Surv(OS.months, Survival.status) ~ HLA_ABCHscore + CD8countpermm2,
               data = HLA_ABCHighonlycuts)
# extract p value
HLA_ABCcd8pvalue <- surv_pvalue(fit)$pval.txt
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  risk.table = FALSE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,60), 
  xlab = ("Months"),# present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 20,    # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), font.x = 25, font.y = 25, font.legend = 20, font.caption = 25,  font.tickslab = 20, pval.size = 10, # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = F, # show bars instead of names in text annotations
  # in legend of risk table
  risk.table.height = 0.3, surv.median.line = "hv",
  legend = "top"
)  + 
  guides(colour = guide_legend(nrow = 2))


#save the image file
dev.copy(png,filename="HLA_ABC_High_Only_CD8OSKm.png");

while (dev.cur()>1) dev.off()
#Does the HLA-ABC H score high only 

library(finalfit)
dependent_os  <- "Surv(OS.months, Survival.status)"
COXphoptimalcutoffsClinMerge$CD8countpermm2
explanatory   <- c(    "Age",
                       "Sex",
                       "pT",
                       "pN",
                       "R0.American.",
                       "MandardScore",
                       "LymphaticInvasion",
                       "VascularInvasion",
                       "HLA_ABCHscore",
                       "CD8countpermm2")

COXphoptimalcutoffsClinMerge %>% 
  finalfit(dependent_os, explanatory)

COXphoptimalcutoffsClinMerge %>% 
  finalfit(dependent_os, explanatory, add_dependent_label = FALSE) %>% 
  rename("Overall survival" = label) %>% 
  rename(" " = levels) %>% 
  rename("  " = all)

COXphoptimalcutoffsClinMerge %>% 
  hr_plot(dependent_os, explanatory)


COXphoptimalcutoffsClinMerge %>%
  coxphmulti(dependent_os, explanatory)
#Assess the HLA-High only + Cd8+ T cell Variables are significant interacting variables
test <- HLA_ABCHighonlyincluded[HLA_ABCHighonlyincluded$HLA_ABC_CD8 !="low:low",]
dependent_os  <- "Surv(OS.months, Survival.status)"
HLA_ABCHighonlyincluded$VascularInvasion
explanatory   <- c(    "Age",
                       "Sex",
                       "pT",
                       "pN",
                       "R0.American.",
                       "MandardScore",
                       "LymphaticInvasion",
                       "VascularInvasion",
                       "HLA_ABC_high_CD8")

HLA_ABCHighonlyincluded %>% 
  finalfit(dependent_os, explanatory)

tab <- HLA_ABCHighonlyincluded %>% 
  finalfit(dependent_os, explanatory, add_dependent_label = FALSE) %>% 
  knitr::kable(row.names=TRUE)
library(rmarkdown)
df_print(tab)
HLA_ABCHighonlyincluded %>% 
  hr_plot(dependent_os, explanatory)


HLA_ABCHighonlyincluded %>%
  coxphmulti(dependent_os, explanatory)




#HLA-ABC CD4 OS survival

fit <- survfit(Surv(OS.months, Survival.status) ~ HLA_ABCHscore + CD4countpermm2,
               data = COXphoptimalcutoffs)
# extract p value
HLA_ABCcd4pvalue <- surv_pvalue(fit)$pval.txt
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,180), 
  xlab = ("Months"),# present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 20,    # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = F, # show bars instead of names in text annotations
  # in legend of risk table
  risk.table.height = 0.3,
  legend = "top"
)  + 
  guides(colour = guide_legend(nrow = 2))


#save the image file
dev.copy(png,filename="HLA_ABCCD4OSKm.png");
while (dev.cur()>1) dev.off()

#HLA-ABC FOXP3 OS survival

fit <- survfit(Surv(OS.months, Survival.status) ~ HLA_ABCHscore + FOXP3countpermm2,
               data = COXphoptimalcutoffs)
# extract p value
HLA_ABCFOXP3pvalue <- surv_pvalue(fit)$pval.txt
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,180), 
  xlab = ("Months"),# present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 20,    # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = F, # show bars instead of names in text annotations
  # in legend of risk table
  risk.table.height = 0.3,
  legend = "top"
)  + 
  guides(colour = guide_legend(nrow = 2))


#save the image file
dev.copy(png,filename="HLA_ABCFOXP3OSKm.png");
while (dev.cur()>1) dev.off()

#HLA-E CD3 OS survival

fit <- survfit(Surv(OS.months, Survival.status) ~ HLA_EHscore + CD3countpermm2,
               data = COXphoptimalcutoffs)
# extract p value
HLA_ECD3pvalue <- surv_pvalue(fit)$pval.txt
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,180), 
  xlab = ("Months"),# present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 20,    # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = F, # show bars instead of names in text annotations
  # in legend of risk table
  risk.table.height = 0.3,
  legend = "top"
)  + 
  guides(colour = guide_legend(nrow = 2))


#save the image file
dev.copy(png,filename="HLA_ECD3OSKm.png");
while (dev.cur()>1) dev.off()

#HLA-E CD8 OS survival

fit <- survfit(Surv(OS.months, Survival.status) ~ HLA_EHscore + CD8countpermm2,
               data = COXphoptimalcutoffs)
# extract p value
HLA_ECD8pvalue <- surv_pvalue(fit)$pval.txt
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,180), 
  xlab = ("Months"),# present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 20,    # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = F, # show bars instead of names in text annotations
  # in legend of risk table
  risk.table.height = 0.3,
  legend = "top"
)  + 
  guides(colour = guide_legend(nrow = 2))


#save the image file
dev.copy(png,filename="HLA_ECD8OSKm.png");
while (dev.cur()>1) dev.off()

#HLA-E CD4 OS survival

fit <- survfit(Surv(OS.months, Survival.status) ~ HLA_EHscore + CD4countpermm2,
               data = COXphoptimalcutoffs)
# extract p value
HLA_ECD4pvalue <- surv_pvalue(fit)$pval.txt
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,180), 
  xlab = ("Months"),# present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 20,    # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = F, # show bars instead of names in text annotations
  # in legend of risk table
  risk.table.height = 0.3,
  legend = "top"
)  + 
  guides(colour = guide_legend(nrow = 2))


#save the image file
dev.copy(png,filename="HLA_ECD4OSKm.png");
while (dev.cur()>1) dev.off()

#HLA-E FOXP3 OS survival

fit <- survfit(Surv(OS.months, Survival.status) ~ HLA_EHscore + FOXp3countpermm2,
               data = COXphoptimalcutoffs)
# extract p value
HLA_EFOXp3pvalue <- surv_pvalue(fit)$pval.txt
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,180), 
  xlab = ("Months"),# present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 20,    # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = F, # show bars instead of names in text annotations
  # in legend of risk table
  risk.table.height = 0.3,
  legend = "top"
)  + 
  guides(colour = guide_legend(nrow = 2))


#save the image file
dev.copy(png,filename="HLA_EFOXP3OSKm.png");
while (dev.cur()>1) dev.off()

#HLA-C2 CD3 OS survival
COXphoptimalcutoffs$HLA_C2Hscore
fit <- survfit(Surv(OS.months, Survival.status) ~ HLA_C2Hscore + CD3countpermm2,
               data = COXphoptimalcutoffs)
# extract p value
HLA_C2CD3pvalue <- surv_pvalue(fit)$pval.txt
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,180), 
  xlab = ("Months"),# present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 20,    # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = F, # show bars instead of names in text annotations
  # in legend of risk table
  risk.table.height = 0.3,
  legend = "top"
)  + 
  guides(colour = guide_legend(nrow = 2))


#save the image file
dev.copy(png,filename="HLA_C2CD3OSKm.png");
while (dev.cur()>1) dev.off()

#HLA-C2 CD8 OS survival
COXphoptimalcutoffs$HLA_C2Hscore
fit <- survfit(Surv(OS.months, Survival.status) ~ HLA_C2Hscore + CD8countpermm2,
               data = COXphoptimalcutoffs)
# extract p value
HLA_C2CD8pvalue <- surv_pvalue(fit)$pval.txt
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,180), 
  xlab = ("Months"),# present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 20,    # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = F, # show bars instead of names in text annotations
  # in legend of risk table
  risk.table.height = 0.3,
  legend = "top"
)  + 
  guides(colour = guide_legend(nrow = 2))


#save the image file
dev.copy(png,filename="HLA_C2CD8OSKm.png");
while (dev.cur()>1) dev.off()

#HLA-C2 CD4 OS survival
COXphoptimalcutoffs$HLA_C2Hscore
fit <- survfit(Surv(OS.months, Survival.status) ~ HLA_C2Hscore + CD4countpermm2,
               data = COXphoptimalcutoffs)
# extract p value
HLA_C2CD4pvalue <- surv_pvalue(fit)$pval.txt
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,180), 
  xlab = ("Months"),# present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 20,    # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = F, # show bars instead of names in text annotations
  # in legend of risk table
  risk.table.height = 0.3,
  legend = "top"
)  + 
  guides(colour = guide_legend(nrow = 2))


#save the image file
dev.copy(png,filename="HLA_C2CD4OSKm.png");
while (dev.cur()>1) dev.off()

#HLA-C2 CD3 OS survival
COXphoptimalcutoffs$HLA_C2Hscore
fit <- survfit(Surv(OS.months, Survival.status) ~ HLA_C2Hscore + CD3countpermm2,
               data = COXphoptimalcutoffs)
# extract p value
HLA_C2CD3pvalue <- surv_pvalue(fit)$pval.txt
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,180), 
  xlab = ("Months"),# present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 20,    # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = F, # show bars instead of names in text annotations
  # in legend of risk table
  risk.table.height = 0.3,
  legend = "top"
)  + 
  guides(colour = guide_legend(nrow = 2))


#save the image file
dev.copy(png,filename="HLA_C2CD3OSKm.png");
while (dev.cur()>1) dev.off()

#heatmap all
#readback the csvs
Allmarkermerged <- read.csv(file = "C:/Users/wp1g19/Documents/TMA analysis/R full analysis wd/AllmarkersmergedQuantileCut.csv", row.names = 1)

library(ALL)
library(pheatmap)
library(RColorBrewer)
library(tidyverse)
breaksList = seq(-3, 3, by = 0.08)
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

omitNAAllmarkerMergedCD <- Allmarkermerged[c(1:4)]
Heat1 <- mutate_all(omitNAAllmarkerMergedCD, function(x) as.numeric(as.character(x)))
Heat1 <- na.omit(Heat1)
Heat2 <- Heat1
Heat2 <- apply(t(Heat2), 1, cal_z_score)
Heat2 <- t(Heat2)
out <- pheatmap(Heat2,
                main = "OAC immune abundances by IHC (69 samples)",
                clustering_distance_rows = "euclidean",
                clustering_method = "ward.D2", 
                clustering_distance_cols = "euclidean",
                cluster_rows = TRUE,
                cluster_cols = TRUE,
                show_colnames = FALSE, cutree_cols = 4,
                breaks = breaksList) 

dev.off()
clusters <- cutree(out$tree_col, k=4)[out$tree_col[["order"]]]
annot_col <- data.frame(col.names = names(clusters),
                        cluster = as.factor(clusters))
annot_col$col.names <- NULL

pheatmap(Heat2,scale = "none",
         main = "OAC immune abundances by IHC (69 samples)",
         clustering_distance_rows = "euclidean",
         clustering_method = "ward.D2", 
         clustering_distance_cols = "euclidean",
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         cutree_cols = 4,
         show_colnames = FALSE,annotation_col = annot_col,
         breaks = breaksList, fontsize_row = 19, fontsize = 18, annotation_legend = TRUE, legend = FALSE)
dev.off()
listA <- c("HLA_ABCHscore", "HLA_C2Hscore", "HLA_EHscore", "TAP1Hscore","CSDE1Hscore")
Allmarkermerged$CSDE1Hscore
AllmarkerMergedCD <- Allmarkermerged[colnames(Allmarkermerged) %in% listA]
AllmarkerMergedCD <- na.omit(AllmarkerMergedCD)
Heat1 <- mutate_all(AllmarkerMergedCD, function(x) as.numeric(as.character(x)))
Heat2 <- Heat1
Heat2 <- apply(t(Heat2), 1, cal_z_score)
Heat2 <- t(Heat2)
out <- pheatmap(Heat2,
                main = "OAC immune abundances by IHC (94 samples)",
                clustering_distance_rows = "euclidean",
                clustering_method = "ward.D2", 
                clustering_distance_cols = "euclidean",
                cluster_rows = TRUE,
                cluster_cols = TRUE,
                show_colnames = FALSE, cutree_cols = 5,
                breaks = breaksList) 

dev.off()
clusters <- cutree(out$tree_col, k=5)[out$tree_col[["order"]]]
annot_col <- data.frame(col.names = names(clusters),
                        cluster = as.factor(clusters))
annot_col$col.names <- NULL
annot_col$cluster <- gsub("5","6", annot_col$cluster)
annot_col$cluster <- gsub("4","7", annot_col$cluster)
annot_col$cluster <- gsub("3","8", annot_col$cluster)
annot_col$cluster <- gsub("1","9", annot_col$cluster)
annot_col$cluster <- gsub("2","10", annot_col$cluster)

annot_col$cluster <- gsub("6","1", annot_col$cluster)
annot_col$cluster <- gsub("7","2", annot_col$cluster)
annot_col$cluster <- gsub("8","3", annot_col$cluster)
annot_col$cluster <- gsub("9","4", annot_col$cluster)
annot_col$cluster <- gsub("10","5", annot_col$cluster)

pheatmap(Heat2,scale = "none",
         main = "OAC APM H score by IHC (94 samples)",
         clustering_distance_rows = "euclidean",
         clustering_method = "ward.D2", 
         clustering_distance_cols = "euclidean",
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         cutree_cols = 5,
         show_colnames = FALSE,annotation_col = annot_col,
         breaks = breaksList, fontsize_row = 19, fontsize = 18, annotation_legend = TRUE, legend = FALSE, border_color = NA)
dev.off()

#merge cluster and TMA data

percentageclusterTMA <- merge(Allmarkermerged, annot_col, by=0)
write.csv(percentageclusterTMA, file = "C:/Users/wp1g19/Documents/TMA analysis/R full analysis wd/images for thesis/Quantile_cut/HscoreAPMclusterTMAmerge.csv")

#Quantile APM expression merge
library(tidyr)
CSDE1naomit <- Allmarkermerged
CSDE1naomit <- CSDE1naomit %>% drop_na(HLA_Epositive.)
quantile(CSDE1naomit$HLA_Epositive.)
CSDE1upperquantile <- CSDE1naomit[CSDE1naomit$HLA_Epositive. > 0.77430000,]
CSDE1lowerquantile <- CSDE1naomit[CSDE1naomit$HLA_Epositive. < 0.01273333,]


#ABC hgh only survival
highonly <- COXphoptimalcutoffsClinMerge


COXphoptimalcutoffsClinMerge$CD8countpermm2
fit <- survfit(Surv(OS.months, Survival.status) ~ HLA_ABCHscore + CD8countpermm2,
               data = COXphoptimalcutoffsClinMerge)
# extract p value

ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,180), 
  xlab = ("Months"),# present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 20,    # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = F, # show bars instead of names in text annotations
  # in legend of risk table
  risk.table.height = 0.3,
  legend = "top"
)  + 
  guides(colour = guide_legend(nrow = 2))
