cormat <- round(cor(as.matrix(omitNAAllmarkerMerged)),2)
cormat <- as.data.frame(cormat)
cormat2 <- cormat[c("CD3countpermm2", "CD4countpermm2", "CD8countpermm2", "FOXP3countpermm2", "HLA_ABCHscore", "HLA_EHscore", "HLA_C2Hscore", "TAP1Hscore"), c("CD3countpermm2", "CD4countpermm2", "CD8countpermm2", "FOXP3countpermm2", "HLA_ABCHscore","HLA_EHscore", "HLA_C2Hscore", "TAP1Hscore" )]
cormat2 <- as.matrix(cormat2)
library(reshape2)
melted_cormat <- melt(cormat2)
head(melted_cormat)
library(ggplot2)
ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()
# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat2){
  cormat2[upper.tri(cormat2)] <- NA
  return(cormat2)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat2){
  cormat2[lower.tri(cormat2)]<- NA
  return(cormat2)
}
upper_tri <- get_upper_tri(cormat2)
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
ggheatmap + 
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1.8, 0.8),
    legend.position = c(1.8, 0.7),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    legend.direction = "vertical", plot.margin = unit(c(1.5,5,0.8,1), "cm"))
dev.copy(png,filename="CorrelationHeatmapOmitNAQuantileCut.png");