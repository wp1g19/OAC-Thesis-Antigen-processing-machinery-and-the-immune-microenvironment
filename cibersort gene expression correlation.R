#read in cibsort data
cibersort <- read.csv(file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Digital cytometry/TPM rerun/CIBERSORTx_Final_Simple.csv", header = TRUE, row.names = 1)
rownames(cibersort) <- gsub("/",".", rownames(cibersort))
rownames(cibersort) <- gsub("-",".", rownames(cibersort))
  #read in gene expression data
  CombinedTMMdata <- read.csv(file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/OCCAMSTCGATMMMergedBatchCorrectedAPMCSDE1OutliersCut.csv", header = TRUE,row.names = 1)
  colnames(CombinedTMMdata) <- strtrim(colnames(CombinedTMMdata), 12)
  rownames(CombinedTMMdata) <- gsub("-",".", rownames(CombinedTMMdata))
  rownames(CombinedTMMdata) <- gsub("/",".", rownames(CombinedTMMdata))
#Filter gene expression to prognostically significant genes
siggenes <- c("RFX5",
              "CTSS",
              "CD1D",
              "MR1",
              "SPPL2A",
              "CD74",
              "CIITA",
              "PSMB10",
              "LGMN",
              "CTSL",
              "TAPBPL",
              "CALR",
              "ERAP1",
              "ERAP2",
              "HLA.E",
              "HLA.DPA1",
              "PSMB9",
              "HLA.B",
              "HLA.DRA",
              "PSMB8",
              "HLA.DRB5",
              "HLA.DQA1",
              "HLA.DRB1",
              "TAPBP",
              "HLA.G",
              "HLA.A",
              "CSDE1")

CombinedTMMdatafilt <- CombinedTMMdata[rownames(CombinedTMMdata) %in% siggenes,]

CombinedTMMdatafilt <- t(CombinedTMMdatafilt)

mergedata <- merge(CombinedTMMdatafilt, cibersort, by = 0)

#filter to survival data
survivaldata <- read.csv(file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/OCCAMSTCGAMergedclinfinal.csv", header = TRUE,row.names = 1)
rownames(survivaldata) <- gsub("/",".", rownames(survivaldata))
rownames(survivaldata) <- gsub("-",".", rownames(survivaldata))
rownames(mergedata) <- mergedata$Row.names
mergedata$Row.names <- NULL
mergedata <- mergedata[rownames(mergedata) %in% rownames(survivaldata),]
#correlation plot


library(psych)
cor_test_mat <- corr.test(mergedata)$p    # Apply corr.test function
cor_test_mat                         # Print matrix of p-values
cor_mat <- cor(mergedata)           # Correlation matrix of example data
cor_mat     
library("corrplot") 
corrplot(cor_mat,                    # Draw corrplot with p-values
         p.mat = cor_test_mat,
         insig = "p-value")
library(ggcorrplot)
ggcorrplot(cor_mat) 

p1 <- ggcorrplot(cor_mat,                  # Draw ggcorrplot with p-values
                 p.mat = cor_test_mat,show.legend = TRUE,method = "square",insig = "pch", lab = TRUE, sig.level = 0.05, pch = 1,pch.cex = 8,pch.col = "red", lab_size = 3 , type = "lower")
p1 + ggtitle("Correlation of MHC II genes and APM regulators n = 176")+ theme(text = element_text(size = 15)) 

p1 <- ggcorrplot(cor_mat,                  # Draw ggcorrplot with p-values
                 p.mat = cor_test_mat,show.legend = TRUE,method = "square",insig = "pch", lab = TRUE, sig.level = 1, pch = 1,pch.cex = 8,pch.col = "red", lab_size = 3 , type = "lower")
p1 + ggtitle("Correlation of MHC II genes and APM regulators n = 176")+ theme(text = element_text(size = 15)) 

#New correlation plot
library(ggplot2)
cor <- Hmisc::rcorr(mergedata %>% as.matrix())
nm = rownames(cor$r)
m = t(combn(nm, 2))
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
  coord_equal()+ theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + ggtitle("Correlation of MHC genes and APM regulators n = 176")+ theme(text = element_text(size = 15)) 

  pvalues <- as.data.frame(cor$P)       
  
  #run after cibersort heatmap script
# Correlation in clusters
  
#read in cibersort data
FRACTcibersort <- read.csv(file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Digital cytometry/TPM rerun/CIBERSORTx_Final_Simple.csv", header = TRUE, row.names = 1)
rownames(FRACTcibersort) <- gsub("/",".", rownames(FRACTcibersort))
rownames(FRACTcibersort) <- gsub("-",".", rownames(FRACTcibersort))
clustersFract <- annot_col
rownames(clustersFract) <- gsub("/",".", rownames(clustersFract))
rownames(clustersFract) <- gsub("-",".", rownames(clustersFract))
cluster1 <- clustersFract[clustersFract$cluster == 1,]
cluster2 <- clustersFract[clustersFract$cluster == 2,]
cluster3 <- clustersFract[clustersFract$cluster == 3,]
cluster4 <- clustersFract[clustersFract$cluster == 4,]

mergedataFract <- merge(CombinedTMMdatafilt, FRACTcibersort, by = 0)
rownames(mergedataFract) <-  mergedataFract$Row.names
mergedataFract$Row.names <- NULL
mergedataFract <- mergedataFract[rownames(mergedataFract) %in% rownames(survivaldata),]

mergedataFractclust1 <- mergedataFract[rownames(mergedataFract) %in% rownames(cluster1),]
rownames(mergedataFractclust1) <-  mergedataFractclust1$Row.names
mergedataFractclust1$Row.names <- NULL
mergedataFractclust2 <- mergedataFract[rownames(mergedataFract) %in% rownames(cluster2),]
rownames(mergedataFractclust2) <-  mergedataFractclust2$Row.names
mergedataFractclust2$Row.names <- NULL
mergedataFractclust3 <- mergedataFract[rownames(mergedataFract) %in% rownames(cluster3),]
rownames(mergedataFractclust3) <-  mergedataFractclust3$Row.names
mergedataFractclust3$Row.names <- NULL
mergedataFractclust4 <- mergedataFract[rownames(mergedataFract) %in% rownames(cluster4),]
rownames(mergedataFractclust4) <-  mergedataFractclust4$Row.names
mergedataFractclust4$Row.names <- NULL

library(psych)
cor_test_mat <- corr.test(mergedataFractclust1)$p    # Apply corr.test function
cor_test_mat                         # Print matrix of p-values
cor_mat <- cor(mergedataFractclust1)           # Correlation matrix of example data
cor_mat     
library("corrplot") 
corrplot(cor_mat,                    # Draw corrplot with p-values
         p.mat = cor_test_mat,
         insig = "p-value")
library(ggcorrplot)
ggcorrplot(cor_mat) 

#Cluster 1 correlation
library(ggplot2)
library(RcmdrMisc)
library(Hmisc)
cor <- rcorr.adjust(mergedataFractclust1 %>% as.matrix(), type = "pearson")
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
  coord_equal()+ theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + ggtitle("Correlation of MHC genes and APM regulators of cluster 1 (n = 9)")+ theme(text = element_text(size = 15)) 

pvalues <- as.data.frame(cor[["R"]][["P"]])       
 write.csv(pvalues, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Digital cytometry/TPM rerun/cluster1sigcorrel.csv")
#Cluster 2 correlation
library(ggplot2)
cor <- Hmisc::rcorr(mergedataFractclust2 %>% as.matrix())
nm = rownames(cor$r)
m = t(combn(nm, 2))
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
  coord_equal()+ theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + ggtitle("Correlation of MHC genes and APM regulators of cluster 2 (n = 9)")+ theme(text = element_text(size = 15)) 

pvalues <- as.data.frame(cor$P)  
write.csv(pvalues, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Digital cytometry/TPM rerun/cluster2sigcorrel.csv")
#Cluster 2 correlation
library(ggplot2)
cor <- Hmisc::rcorr(mergedataFractclust3 %>% as.matrix())
nm = rownames(cor$r)
m = t(combn(nm, 2))
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
  coord_equal()+ theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + ggtitle("Correlation of MHC genes and APM regulators of cluster 3 (n = 40)")+ theme(text = element_text(size = 15)) 
pvalues <- as.data.frame(cor$P) 
write.csv(pvalues, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Digital cytometry/TPM rerun/cluster3sigcorrel.csv")
#Cluster 2 correlation
library(ggplot2)
cor <- Hmisc::rcorr(mergedataFractclust4 %>% as.matrix())
nm = rownames(cor$r)
m = t(combn(nm, 2))
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
  coord_equal()+ theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + ggtitle("Correlation of MHC genes and APM regulators of cluster 4 (n = 103)")+ theme(text = element_text(size = 15)) 

pvalues <- as.data.frame(cor$P) 
write.csv(pvalues, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Digital cytometry/TPM rerun/cluster4sigcorrel.csv")
#correlation between genes and absolute score
absolutescore <- read.csv("C:/Users/wp1g19/Documents/Revisit bioinformatics/Digital cytometry/TPM rerun/ABSOLUTESCORE.csv", header = TRUE, row.names = 1)
rownames(absolutescore) <- gsub("/",".", rownames(absolutescore))
rownames(absolutescore) <- gsub("-",".", rownames(absolutescore))
mergedataTIL <- merge(CombinedTMMdatafilt, absolutescore, by = 0)
rownames(mergedataTIL) <- mergedataTIL$Row.names
mergedataTIL$Row.names <- NULL

#TILS correlation
library(ggplot2)
cor <- Hmisc::rcorr(mergedataTIL %>% as.matrix())
nm = rownames(cor$r)
m = t(combn(nm, 2))
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
  coord_equal()+ theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + ggtitle("Correlation of MHC genes and ABSOLUTE TILs") + theme(text = element_text(size = 15)) 

