#ECOTYPE/CIBESORT

#READ in cutpoints and format (OS cutpoints already filtered to significance)
OScutpoints <- read.csv(file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/optimalcutpointsOS.csv", header = TRUE, row.names = 1)
OScutpoints$ID <- NULL
CSScutpoints <- read.csv(file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/optimalcutpointsCSS.csv", header = TRUE, row.names = 1)
CSScutpoints$OS.months <- NULL
CSScutpoints$Bioinformatic.CSS_status <- NULL
DFScutpoints <- read.csv(file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/optimalcutpointsDFS.csv", header = TRUE, row.names = 1)
DFScutpoints$Bioinformatic.DFS_months <- NULL
DFScutpoints$Bioinformatic.DFS_status <- NULL

#read in ecotyper carcinoma TIME 
ecotypersCarciomaTIME <- read.csv("C:/Users/wp1g19/Documents/Revisit bioinformatics/Digital cytometry/ecotyper_output ALL/Carcinoma_Ecotypes/Ecotype_Abundance.csv", header = TRUE, row.names = 1)
rownames(ecotypersCarciomaTIME) <- strtrim(rownames(ecotypersCarciomaTIME),12)

#filter ecotyper data to OS,CSS,DFS data 
ecotypersCarciomaTIMEfiltOS <- ecotypersCarciomaTIME[rownames(ecotypersCarciomaTIME) %in% rownames(OScutpoints),]
ecotypersCarciomaTIMEfiltCSS <- ecotypersCarciomaTIME[rownames(ecotypersCarciomaTIME) %in% rownames(CSScutpoints),]
ecotypersCarciomaTIMEfiltDFS <- ecotypersCarciomaTIME[rownames(ecotypersCarciomaTIME) %in% rownames(DFScutpoints),]

#Merge the ecotyper TIME type data with the cutpoints
ecotypersCarciomaTIMEfiltOS$ID <- rownames(ecotypersCarciomaTIMEfiltOS)
OScutpoints$ID <- row.names(OScutpoints)
mergeECO_OS <- merge.data.frame(OScutpoints, ecotypersCarciomaTIMEfiltOS, by = "ID")
ecotypersCarciomaTIMEfiltCSS$ID <- rownames(ecotypersCarciomaTIMEfiltCSS)
CSScutpoints$ID <- row.names(CSScutpoints)
mergeECO_CSS <- merge.data.frame(CSScutpoints, ecotypersCarciomaTIMEfiltCSS, by = "ID")
ecotypersCarciomaTIMEfiltDFS$ID <- rownames(ecotypersCarciomaTIMEfiltDFS)
DFScutpoints$ID <- row.names(DFScutpoints)
mergeECO_DFS <- merge.data.frame(DFScutpoints, ecotypersCarciomaTIMEfiltDFS, by = "ID")

#Remove ID column
rownames(mergeECO_OS) <- mergeECO_OS$ID
rownames(mergeECO_CSS) <- mergeECO_CSS$ID
rownames(mergeECO_DFS)<- mergeECO_DFS$ID
mergeECO_OS$ID <- NULL
mergeECO_CSS$ID <- NULL
mergeECO_DFS$ID <- NULL

# Create data frame with ecotype assignment of maximum abundance of type
ECOassignOS <- as.data.frame(row.names(ecotypersCarciomaTIMEfiltOS))
ECOassignOS$Ecotype <- colnames(ecotypersCarciomaTIMEfiltOS)[apply(ecotypersCarciomaTIMEfiltOS,1,which.max)]
colnames(ECOassignOS) <- c("ID", "Ecotype")
write.csv(ECOassignOS, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Digital cytometry/ECOassignmentOS.csv")

#heatmap the ecotype states
library(dplyr)
ecotypersCarciomaTIMEfiltOS$ID <- NULL
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}
ecotypersCarciomaTIMEfiltOSHEAT <- mutate_all(ecotypersCarciomaTIMEfiltOS, function(x) as.numeric(as.character(x)))
ecotypersCarciomaTIMEfiltOSHEAT <- apply(ecotypersCarciomaTIMEfiltOS, 1, cal_z_score)
#ecotypersCarciomaTIMEfiltOSHEAT <- t(ecotypersCarciomaTIMEfiltOSHEAT)
library(ALL)
library(pheatmap)
library(RColorBrewer)
library(tidyverse)
breaksList = seq(-3, 3, by = 0.08)
annot_col <- ECOassignOS
rownames(annot_col) <- annot_col$ID
annot_col$ID <- NULL
annot_col %>% arrange(annot_col)
ecotypersCarciomaTIMEfiltOSHEAT <- ecotypersCarciomaTIMEfiltOSHEAT[, rownames(annot_col) ]
out <- pheatmap(ecotypersCarciomaTIMEfiltOSHEAT,scale = "none",
                main = "OAC Ecotypes (176 samples)",
                clustering_distance_rows = "euclidean",
                clustering_method = "ward.D2", 
                clustering_distance_cols = "euclidean",
                cluster_rows = FALSE,
                cluster_cols = TRUE,
                show_colnames = FALSE,annotation_col = annot_col,
                breaks = breaksList) 
dev.off()

out <- pheatmap(ecotypersCarciomaTIMEfiltOSHEAT,scale = "none",
                main = "OAC Ecotypes (176 samples)",
                clustering_distance_rows = "euclidean",
                clustering_method = "ward.D2", 
                clustering_distance_cols = "euclidean",
                cluster_rows = FALSE,
                cluster_cols = TRUE,
                show_colnames = FALSE,cutree_cols = 10,annotation_col = annot_col,
                breaks = breaksList) 
