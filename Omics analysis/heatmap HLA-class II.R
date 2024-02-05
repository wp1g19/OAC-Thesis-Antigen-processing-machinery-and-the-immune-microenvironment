batchedTCGAOCCAMS <- read.csv(file ="C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/OCCAMSTCGATMMMergedBatchCorrected.csv", header = TRUE, row.names = 1)
APMgenes <- c(     "HLA-DMA",
                   "HLA-DRA",
                   "HLA-DPA1",
                   "HLA-DRB5",
                   "HLA-DQA1",
                   "HLA-DRB1",
                   "HLA-DOA",
                   "HLA-DQA2"
                   )


batchedTCGAOCCAMSAPM <- batchedTCGAOCCAMS[ rownames(batchedTCGAOCCAMS) %in% APMgenes ,]
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}
batchedTCGAOCCAMSAPM <- apply(batchedTCGAOCCAMSAPM, 1, cal_z_score)
batchedTCGAOCCAMSAPM <- t(batchedTCGAOCCAMSAPM)
library(ALL)
library(pheatmap)
library(RColorBrewer)
breaksList = seq(-3, 3, by = 0.08)
out <- pheatmap(batchedTCGAOCCAMSAPM,scale = "none",
                main = "APM heatmap in OAC cohort (194 samples, TMM)",
                clustering_distance_rows = "euclidean",
                clustering_method = "ward.D2",
                clustering_distance_cols = "euclidean",
                cluster_rows = TRUE,
                cluster_cols = TRUE,
                cutree_cols = 3,
                show_colnames = FALSE,
                breaks = breaksList) 
clusters <- cutree(out$tree_col, k=3)[out$tree_col[["order"]]]
annot_col <- data.frame(col.names = names(clusters),
                        cluster = as.factor(clusters))
annot_col$col.names <- NULL


pheatmap(batchedTCGAOCCAMSAPM,scale = "none",
         main = "APM heatmap in OAC cohort (194 samples, TMM)",
         clustering_distance_rows = "euclidean",
         clustering_method = "ward.D2", 
         clustering_distance_cols = "euclidean",
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         cutree_cols = 3,
         show_colnames = FALSE,annotation_col = annot_col,
         breaks = breaksList)

Clusters <- as.data.frame(clusters)

survdata <- read.csv("C:/Users/wp1g19/OneDrive - University of Southampton/General/TCGAOCCAMS/OCCAMSTCGAMergedclinfinal.csv")

Clusters$ID <- strtrim(row.names(Clusters), 12)
Clusters$ID <- gsub('\\.', '/', Clusters$ID)
survdata$ID <- gsub('\\.', '/', survdata$ID)
clustersSurvmerge <- merge(Clusters, survdata, by.x = "ID")


#survival on clusters

library("tidyverse")
library(tidyr)
library(survival)
library(survminer)
library(RTCGA)
library(forestmodel)
library(finalfit)
library(survivalAnalysis)

clustersSurvmerge$Vital_status <- gsub('yes', '0', clustersSurvmerge$Vital_status)
clustersSurvmerge$Vital_status <- gsub('no', '1', clustersSurvmerge$Vital_status)
clustersSurvmerge$Vital_status <- as.numeric(clustersSurvmerge$Vital_status)
fit <- survfit(Surv(OS.months, Vital_status) ~ clusters,
               data = clustersSurvmerge)

ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,40), 
  xlab = ("Months"),# present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 10,    # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = TRUE, surv.median.line = "hv" # show bars instead of names in text annotations
  # in legend of risk table
)  
