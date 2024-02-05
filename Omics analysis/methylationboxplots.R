library(TCGAbiolinks)
library(tidyverse)

TCGAesca <- GDCquery_clinic("TCGA-ESCA")
setwd("C:/Users/wp1g19/Documents/Revisit bioinformatics/Methlyation/")
##Read files 


data_files <- list.files("C:/Users/wp1g19/Documents/Revisit bioinformatics/Methlyation")  # Identify file names
data_files                                                 
for(i in 1:length(data_files)) {                              # Head of for-loop
  assign(paste0("data", i),                                   # Read and store data frames
         read.csv(paste0("C:/Users/wp1g19/Documents/Revisit bioinformatics/Methlyation/",
                          data_files[i])))
}

df_list <- list(data1,data2,data3,data4, data5,data6, data7, data8,data9,data10)
cpgsmerge <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list)
colnames(cpgsmerge) <- c( "X","sample", "tissue", "type", "CALR", "CIITA", "CSDE1", "ERAP2", "HLA-A", "HLA-B", "HLA-DPA1", "HLA-DPB1", "HLA-G", "RFX5")
cpgsmerge$sample <- strtrim(cpgsmerge$sample,12)
splitdf <- split.data.frame(cpgsmerge, f = cpgsmerge$type)

Tumourcpgs <- splitdf$Tumor
Normalcpgs <- splitdf$Normal
colnames(TCGAesca)[2] ="sample"
TUMOURCPGmerge <- left_join(Tumourcpgs, TCGAesca, by= "sample")
NORMALCPGmerge <- left_join(Normalcpgs, TCGAesca, by= "sample")

OACCPGmerge <- TUMOURCPGmerge[TUMOURCPGmerge$primary_diagnosis == "Adenocarcinoma, NOS",]
NORMALOACCPGmerge <- NORMALCPGmerge[NORMALCPGmerge$primary_diagnosis == "Adenocarcinoma, NOS",]

OAC_Normalcpgsmerged <- rbind(OACCPGmerge, NORMALOACCPGmerge)

write.csv(OAC_Normalcpgsmerged, "C:/Users/wp1g19/Documents/Revisit bioinformatics/Methlyation/OAC_Normalcpgsmerged.csv")
library(hrbrthemes)
library(viridis)
library(ggpubr)
OAC_Normalcpgsmerged %>%
  ggplot( aes(x=type, y=`HLA-A`, fill=type)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("HLA-A methylations Normal vs OAC") +
  xlab("")+ stat_compare_means(method = "t.test")

OAC_Normalcpgsmerged %>%
  ggplot( aes(x=type, y=`HLA-B`, fill=type)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("HLA-B methylations Normal vs OAC") +
  xlab("")+ stat_compare_means(method = "t.test")

OAC_Normalcpgsmerged %>%
  ggplot( aes(x=type, y=`HLA-G`, fill=type)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("HLA-G methylations Normal vs OAC") +
  xlab("")+ stat_compare_means(method = "t.test")

OAC_Normalcpgsmerged %>%
  ggplot( aes(x=type, y=`HLA-DPA1`, fill=type)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("HLA-DPA1 methylations Normal vs OAC") +
  xlab("")+ stat_compare_means(method = "t.test")

OAC_Normalcpgsmerged %>%
  ggplot( aes(x=type, y=`HLA-DPB1`, fill=type)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("HLA-DPB1 methylations Normal vs OAC") +
  xlab("")+ stat_compare_means(method = "t.test")

OAC_Normalcpgsmerged %>%
  ggplot( aes(x=type, y= ERAP2, fill=type)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("ERAP2 methylations Normal vs OAC") +
  xlab("")+ stat_compare_means(method = "t.test")

OAC_Normalcpgsmerged %>%
  ggplot( aes(x=type, y=`CSDE1`, fill=type)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("CSDE1 methylations Normal vs OAC") +
  xlab("")+ stat_compare_means(method = "t.test")

OAC_Normalcpgsmerged %>%
  ggplot( aes(x=type, y=`CIITA`, fill=type)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("CIITA methylations Normal vs OAC") +
  xlab("")+ stat_compare_means(method = "t.test")

OAC_Normalcpgsmerged %>%
  ggplot( aes(x=type, y=`RFX5`, fill=type)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("RFX5 methylations Normal vs OAC") +
  xlab("")+ stat_compare_means(method = "t.test")

OAC_Normalcpgsmerged %>%
  ggplot( aes(x=type, y=`CALR`, fill=type)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("CALR methylations Normal vs OAC") +
  xlab("")+ stat_compare_means(method = "t.test")
OACCPGmerge$ajcc_pathologic_t
#Stageplots
OACCPGmerge %>%
  ggplot( aes(x=ajcc_pathologic_t, y=`CSDE1`, fill=ajcc_pathologic_t)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("CALR methylations Normal vs OAC") +
  xlab("")+ stat_compare_means(method = "anova")


