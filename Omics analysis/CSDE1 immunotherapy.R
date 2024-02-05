library(edgeR)
library(dplyr)
library(tidyr)
#ICIcounts <- read.table(file="C:/Users/wp1g19/Documents/Revisit bioinformatics/CSDE1 in immunotherapy responses/iatlas-ici-hgnc_counts.tsv", header = TRUE, row.names =  1)
#write.csv(ICIcounts, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/CSDE1 in immunotherapy responses/iatlas-ici-hgnc_counts.csv")
ICIcounts <- read.csv(file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/CSDE1 in immunotherapy responses/iatlas-ici-hgnc_counts.csv", header = TRUE, row.names =  1)
ICIcounts <- t(ICIcounts)
dge <- DGEList(counts = ICIcounts)
keep  <- filterByExpr(dge)
table(keep)
dge <- dge[keep, , keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge,method = "TMM", lib.size=NULL)
ICITMMMatrix <- cpm(dge)

pca <- mixOmics::pca(t(ICITMMMatrix), ncomp = 3)
colorlist <- rainbow(12)
library(PLSDAbatch)
Scatter_Density(pca)
library(edgeR)
dev.off()
batch <- c(rep("Choueiri_2016",16),
           rep("Gide_2019", 91),
           rep("HugoLo_2016", 27),
           rep("IMmotion150", 263),
           rep("IMVigor210", 348),
           rep("Kim_2018", 78),
           rep("Liu_2019", 122),
           rep("Miao_2018", 17),
           rep("Prins_2019", 30),
           rep("Riaz_2017", 107),
           rep("VanAllen_2015", 42),
           rep("Zhao_2019", 34))
Scatter_Density(pca, batch = batch, color.set = colorlist)
dev.off()
library(sva)
ICITMMMatrixBatch <- ComBat(
  ICITMMMatrix,
  batch ,
  mod = NULL,
  par.prior = TRUE,
  prior.plots = FALSE,
  mean.only = FALSE,
  ref.batch = NULL,
  BPPARAM = bpparam("SerialParam")
)
colsnam <- colnames(ICITMMMatrixBatch)
pca <- mixOmics::pca(t(ICITMMMatrixBatch), ncomp = 3)
Scatter_Density(pca, batch = batch, color.set = colorlist)
dev.off()
sampledata <- read.csv("C:/Users/wp1g19/Documents/Revisit bioinformatics/CSDE1 in immunotherapy responses/iatlas-ici-sample_info.csv", header = TRUE, row.names = 1)

genes <- c("CSDE1", "HLA.A", "HLA.B", "HLA.C", "HLA.E", "CD274", "PDCD1")

ICITMMMatrixBatchFilt <- as.data.frame(ICITMMMatrixBatch[rownames(ICITMMMatrixBatch) %in% genes ,])
colnames(ICITMMMatrixBatchFilt) <- c("CSDE1", "HLA.A", "HLA.B", "HLA.C", "HLA.E", "PD-L1", "PD-1")
ICITMMMatrixBatchFilt <- as.data.frame(t(ICITMMMatrixBatchFilt))
ICITMMMatrixBatchFilt$RunID <- colnames(ICITMMMatrixBatch)
sampledata$RunID <- rownames(sampledata)
mergedata <- left_join(sampledata, ICITMMMatrixBatchFilt, "RunID")

library(tidyverse)
library(hrbrthemes)
library(viridis)


mergedata %>%
  ggplot( aes(x=Response, y=CSDE1, fill=Response)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("A boxplot with jitter") +
  xlab("")


write.csv(mergedata,file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/CSDE1 in immunotherapy responses/CSDE1Immunotherapymerged.csv")

response <- c("Partial Response", "Complete Response")
responders <- as.data.frame(mergedata[mergedata$Response %in% response ,])
Noresponse <- c("Progressive Disease", "Stable Disease")
Nonresponders <- as.data.frame(mergedata[mergedata$Response %in% Noresponse ,])

write.csv(responders,file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/CSDE1 in immunotherapy responses/CSDE1respondersImmunotherapymerged.csv")
write.csv(Nonresponders,file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/CSDE1 in immunotherapy responses/CSDE1NonrespondersImmunotherapymerged.csv")



# linear trend + confidence interval
library(ggplot2)
library(ggpmisc)
p1 <- ggplot(mergedata, aes(x=mergedata$CSDE1, y=mergedata$PDCD1)) +
  geom_point() +
  geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
  theme_ipsum()+
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "f", "p", "n"))) +
  geom_point()
p1

p2 <- ggplot(mergedata, aes(x=mergedata$CSDE1, y=mergedata$CD274)) +
  geom_point() +
  geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
  theme_ipsum() +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "f", "p", "n"))) +
  geom_point()
p2
