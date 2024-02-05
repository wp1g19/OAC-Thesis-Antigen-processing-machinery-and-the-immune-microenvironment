#scree plot PCA immune
Allmarkermerged <- read.csv(file = "C:/Users/wp1g19/Documents/TMA analysis/R full analysis wd/AllmarkersmergedQuantileCut.csv", row.names = 1)
immune <- Allmarkermerged[c(1:4)]
immune <- na.omit(immune)
#perform PCA
results <- prcomp(immune, scale = TRUE)
#calculate total variance explained by each principal component
var_explained = results$sdev^2 / sum(results$sdev^2)

#create scree plot
library(ggplot2)

qplot(c(1:4), var_explained) + 
  geom_line() + 
  xlab("Principal Component") + 
  ylab("Variance Explained") +
  ggtitle("Scree Plot") +
  ylim(0, 1)
print(var_explained)
# Hscores
Hscore <- Allmarkermerged[c("HLA_ABCHscore", "HLA_C2Hscore", "HLA_EHscore", "TAP1Hscore", "CSDE1Hscore")]
Hscore <- na.omit(Hscore)
#perform PCA
results <- prcomp(Hscore, scale = TRUE)
#calculate total variance explained by each principal component
var_explained = results$sdev^2 / sum(results$sdev^2)

#create scree plot
library(ggplot2)

qplot(c(1:5), var_explained) + 
  geom_line() + 
  xlab("Principal Component") + 
  ylab("Variance Explained") +
  ggtitle("Scree Plot") +
  ylim(0, 1)
print(var_explained)

# Hscores

perc <- Allmarkermerged[c("HLA_ABCpercentagescore", "HLA_C2positive.", "HLA_Epositive.", "TAP1positive.", "CSDE1positive.")]
perc <- na.omit(perc)
#perform PCA
results <- prcomp(perc, scale = TRUE)
#calculate total variance explained by each principal component
var_explained = results$sdev^2 / sum(results$sdev^2)

#create scree plot
library(ggplot2)

qplot(c(1:5), var_explained) + 
  geom_line() + 
  xlab("Principal Component") + 
  ylab("Variance Explained") +
  ggtitle("Scree Plot") +
  ylim(0, 1)
print(var_explained)
library(tidyverse)
library(magrittr)
library(cluster)
library(cluster.datasets)
library(cowplot)
library(NbClust)
library(clValid)
library(ggfortify)
library(clustree)
library(dendextend)
library(factoextra)
library(FactoMineR)
library(corrplot)
library(GGally)
library(ggiraphExtra)
library(knitr)
library(kableExtra)
res.pca <- PCA(Hscore,  graph = FALSE)
# Visualize eigenvalues/variances

fviz_nbclust(Hscore, kmeans, method = "wss", k.max = 12) + theme_minimal() + ggtitle("the Elbow Method")
fviz_nbclust(immune, kmeans, method = "wss", k.max = 12) + theme_minimal() + ggtitle("the Elbow Method")
fviz_nbclust(Hscore, kmeans, method = "wss", k.max = 12) + theme_minimal() + ggtitle("the Elbow Method")

hscore2 <- as.matrix(Hscore)
hscore3 <- clValid(hscore2, nClust = 5:24, 
                  clMethods = c("hierarchical","kmeans","pam"), validation = "internal")
# Summary
summary(hscore3) %>% kable() %>% kable_styling()

immune1 <- as.matrix(immune)
immune2 <- clValid(immune1, nClust = 4:24, 
                   clMethods = c("hierarchical","kmeans","pam"), validation = "internal")
# Summary
summary(immune2) %>% kable() %>% kable_styling()

perc1 <- as.matrix(perc)
perc2 <- clValid(perc1, nClust = 5:24, 
                   clMethods = c("hierarchical","kmeans","pam"), validation = "internal")
# Summary
summary(immune2) %>% kable() %>% kable_styling()

