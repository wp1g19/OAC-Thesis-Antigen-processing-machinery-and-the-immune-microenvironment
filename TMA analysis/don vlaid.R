deconvalid <- read.csv(file = "C:/Users/wp1g19/Documents/TMA analysis/R full analysis wd/Deconvolution validation.csv", row.names = 1)

#P adjusted heatmap correlation

library(ggplot2)
library(RcmdrMisc)
library(Hmisc)
cor <- rcorr.adjust(deconvalid %>% as.matrix(), type = "pearson")
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

#TCGA
#read in cibsort data
library(dplyr)
cibersort <- read.csv(file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Digital cytometry/TPM rerun/CIBERSORTx_Final_Simple.csv", header = TRUE, row.names = 1)
rownames(cibersort) <- gsub("/",".", rownames(cibersort))
rownames(cibersort) <- gsub("-",".", rownames(cibersort))

escabiospecimen <- read.csv(file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Digital cytometry/TCGA-ESCA biospecimen.csv", header = TRUE)
escabiospecimen$ID <- gsub("/",".", escabiospecimen$ID)
escabiospecimen$ID <- gsub("-",".", escabiospecimen$ID)
escabiospecimen2 <- group_by(escabiospecimen, ID) %>% summarize(m = mean(percent_lymphocyte_infiltration))
colnames(escabiospecimen2)  <- c("ID", "percent_lymphocyte_infiltration")
cibersort2 <- cibersort[c(3,4,5)]
cibersort2 <- cibersort2 %>%
  mutate(TILsabsolutescore = rowSums(.))
print(cibersort2)

cibersort3 <- cibersort2[ rownames(cibersort2) %in% escabiospecimen2$ID ,]
escabiospecimen2 <- escabiospecimen2[ escabiospecimen2$ID %in% rownames(cibersort3),]
cibersort3$ID <- rownames(cibersort3)

cibersortbiospecimenmerge <- merge(cibersort3, escabiospecimen2, by = "ID")
