#methylation expression correlation
methylationdata <- read.csv("C:/Users/wp1g19/Documents/Revisit bioinformatics/Methlyation/OACcpgsmerged.csv", header = TRUE)

expressiondata <- read.csv("C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/OCCAMSTCGATMMMergedBatchCorrected.csv", header = TRUE, row.names = 1)


expressiondata <- t(expressiondata)
expressiondata <- as.data.frame(expressiondata)
expressiondata$sample <- rownames(expressiondata)
expressiondata$sample <- strtrim(expressiondata$sample, 12)
methylationdata$sample <- gsub("-", ".", methylationdata$sample)
rownames(expressiondata) <- expressiondata$sample
rownames(methylationdata) <- methylationdata$sample
COLS <- c("sample",
"CALR",                                              
"CIITA",
"CSDE1", 
"ERAP2",
"HLA.A",                                              
"HLA.B",
"HLA.DPA1",                                            
"HLA.DPB1",
"HLA.G",                                               
"RFX5")             
colnames(expressiondata) <- gsub("-", ".", colnames(expressiondata))

expressiondata <- expressiondata[colnames(expressiondata) %in% COLS]
methylationdata <- methylationdata[colnames(methylationdata) %in% COLS]
expressiondata <- expressiondata[expressiondata$sample %in% rownames(methylationdata),]
methylationdata <- methylationdata[methylationdata$sample %in% rownames(expressiondata),]


Methexprmerge <- as.data.frame(merge(methylationdata, expressiondata, by ="sample"))


colsname <- c("sample"  ,
"CALR.Methylation.B.value",
"CIITA.Methylation.B.value",
"CSDE1.Methylation.B.value",
"ERAP2.Methylation.B.value" ,  
"HLA.A.Methylation.B.value",
"HLA.B.Methylation.B.value",
"HLA.DPA1.Methylation.B.value",
"HLA.DPB1.Methylation.B.value",
"HLA.G.Methylation.B.value",   
"RFX5.Methylation.B.value",
"CSDE1.Expression.TMM",
"RFX5.Expression.TMM",
"ERAP2.Expression.TMM",
"CALR.Expression.TMM",
"CIITA.Expression.TMM", 
"HLA.G.Expression.TMM",
"HLA.A.Expression.TMM" ,
"HLA.DPB1.Expression.TMM", 
"HLA.DPA1.Expression.TMM" ,                 
"HLA.B.Expression.TMM" )

colnames(Methexprmerge) <- colsname
library(ggplot2)
library(hrbrthemes)
ggplot(Methexprmerge, aes(x=CSDE1.Methylation.B.value, y=`CSDE1.Expression.TMM`)) +
  geom_point() +
  geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
  stat_cor(method = "pearson", label.x = 0, label.y = 1)+ theme(text = element_text(size = 24)) 

ggplot(Methexprmerge, aes(x=CALR.Methylation.B.value, y=CALR.Expression.TMM)) +
  geom_point() +
  geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
  theme_ipsum()+
  stat_cor(method = "pearson", label.Methylation.B.value = 0, label.y = 1)

ggplot(Methexprmerge, aes(x=CIITA.Methylation.B.value, y=CIITA.Expression.TMM)) +
  geom_point() +
  geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
  stat_cor(method = "pearson", label.x = 0, label.y = 1)+ theme(text = element_text(size = 24)) 

ggplot(Methexprmerge, aes(x=ERAP2.Methylation.B.value, y=ERAP2.Expression.TMM)) +
  geom_point() +
  geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
  theme_ipsum()+
  stat_cor(method = "pearson", label.x = 0, label.y = 1)

ggplot(Methexprmerge, aes(x=HLA.A.Methylation.B.value, y=HLA.A.Expression.TMM)) +
  geom_point() +
  geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
  theme_ipsum()+
  stat_cor(method = "pearson", label.x = 0, label.y = 1)

ggplot(Methexprmerge, aes(x=HLA.B.Methylation.B.value, y=HLA.B.Expression.TMM)) +
  geom_point() +
  geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
  stat_cor(method = "pearson", label.x = 0, label.y = 1)+ theme(text = element_text(size = 24)) 

ggplot(Methexprmerge, aes(x=HLA.DPA1.Methylation.B.value, y=HLA.DPA1.Expression.TMM)) +
  geom_point() +
  geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
  theme_ipsum()+
  stat_cor(method = "pearson", label.x = 0, label.y = 1)

ggplot(Methexprmerge, aes(x=HLA.DPB1.Methylation.B.value, y=HLA.DPB1.Expression.TMM)) +
  geom_point() +
  geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
  theme_ipsum()+
  stat_cor(method = "pearson", label.x = 0, label.y = 1)

ggplot(Methexprmerge, aes(x=HLA.G.Methylation.B.value, y=HLA.G.Expression.TMM)) +
  geom_point() +
  geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
  stat_cor(method = "pearson", label.x = 0, label.y = 1)+ theme(text = element_text(size = 24)) 

ggplot(Methexprmerge, aes(x=RFX5.Methylation.B.value, y=RFX5.Expression.TMM)) +
  geom_point() +
  geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
  theme_ipsum()+
  stat_cor(method = "pearson", label.x = 0, label.y = 1)
