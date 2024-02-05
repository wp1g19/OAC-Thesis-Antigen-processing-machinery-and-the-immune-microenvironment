library(gplots)
library(Heatplus)
library(dendextend)
library(pheatmap)

#-----------------------
# combined deconvolution
optimalcuts <- read.csv(file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/CombinedOSoptimalcut.csv", header=T, stringsAsFactors=F, row.names = 1)

cibersortabs <- read.csv(file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/AbsoluteCibersortRevisit.csv", header=T, stringsAsFactors=F, row.names = 1)
cibersortfractional <- read.csv(file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/CibersortTCGAOCCAMSrevisitFractional.csv", header=T, stringsAsFactors=F, row.names = 1)
rownames(cibersortabs) <- gsub("/", ".", rownames(cibersortabs))
rownames(cibersortabs) <- gsub("-", ".", rownames(cibersortabs))
rownames(cibersortfractional) <- gsub("/", ".", rownames(cibersortfractional))
rownames(cibersortfractional) <- gsub("-", ".", rownames(cibersortfractional))
cm <- as.matrix(cibersortabs) # column 6 is all 0's

exprs <- cm

#calculate z scores
calc_z_score <- function(X,
                         mean,
                         sd) {
  X_scale <- matrix(0, nrow(X), ncol(X),
                    dimnames = list(rownames(X), colnames(X))
  )
  
  if (missing(mean) & missing(sd)) {
    mean_X <- colMeans(X, na.rm = TRUE)
    sd_X <- matrixStats::colSds(as.matrix(X), na.rm = TRUE)
    X_scale <- sweep(X, 2, mean_X, FUN = "-")
    X_scale <- sweep(X_scale, 2, sd_X, FUN = "/")
  } else {
    mean <- mean[na.omit(match(colnames(X), names(mean)))]
    sd <- sd[na.omit(match(colnames(X), names(sd)))]
    
    X <- X[, na.omit(match(names(sd), colnames(X)))]
    X_scale <- sweep(X, 2, mean, FUN = "-")
    X_scale <- sweep(X_scale, 2, sd, FUN = "/")
  }
  return(as.matrix(X_scale))
}


matOAC <- calc_z_score(t(exprs))
library(pheatmap)
pheatmap(matOAC, show_colnames = FALSE, cluster_rows = FALSE, cluster_cols = FALSE)

pheatmap(matOAC, show_colnames = FALSE, cluster_rows = TRUE, cluster_cols = TRUE, scale = "none")

pheatmap(matOAC, show_rownames = TRUE, show_colnames = FALSE, cluster_rows = TRUE, cluster_cols = TRUE, scale = "none", color = colorRampPalette(c("navy", "white", "red"))(50))

anno <- read.csv("C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/absolutecibersortanno.csv", header=T, stringsAsFactors=F, row.names = 1)
my_sample_col <- data.frame(anno)

pheatmap(matOAC,main = "Combined deconvolution n= 179",cutree_cols = 4, clustering_method = "complete", show_rownames = TRUE, show_colnames = FALSE, cluster_rows = TRUE, cluster_cols = TRUE, scale = "none",annotation_col = my_sample_col, color = colorRampPalette(c("navy", "white", "red"))(100))

# apply high low expression groups to cibersort in correlation, boxplot and heatmaps

#start with heatmaps
optimalcuts2 <- optimalcuts
optimalcuts2$OS.months <-NULL
optimalcuts2$Vital_status <- NULL
my_sample_col2 <- data.frame(optimalcuts2)
my_sample_col2 <- my_sample_col2[rownames(my_sample_col2) %in%  rownames(cibersortabs),]
genes <- c("CSDE1", "HLA_A", "ERAP2", "MR1", "LGMN")
my_sample_col3 <- as.data.frame(my_sample_col2[,colnames(my_sample_col2) %in%  genes])



pheatmap(matOAC,main = "Combined deconvolution n= 179",cutree_cols = 4, clustering_method = "complete", show_rownames = TRUE, show_colnames = FALSE, cluster_rows = TRUE, cluster_cols = TRUE, scale = "none",annotation_col = my_sample_col3, color = colorRampPalette(c("navy", "white", "red"))(100))

#Boxplot

boxdata <- merge(my_sample_col2, cibersortabs, by = 0)

X<-split(boxdata, boxdata$CSDE1)

boxdata2 <- X$high
boxdata3 <- X$low

cibersimple <- read.csv("C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/SimpleCibersortsort.csv", header=T, stringsAsFactors=F, row.names = 1)
rownames(cibersimple) <- gsub("/", ".", rownames(cibersimple))
rownames(cibersimple) <- gsub("-", ".", rownames(cibersimple))

Cibersimplehigh <- as.data.frame(cibersimple[rownames(cibersimple) %in%  boxdata2$Row.names,])
Cibersimplelow <- as.data.frame(cibersimple[rownames(cibersimple) %in%  boxdata3$Row.names,])

write.csv(Cibersimplehigh, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/SimpleCibersortsortcsde1high.csv")
write.csv(Cibersimplelow, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/SimpleCibersortsortcsde1low.csv")


write.csv(UpperquantileCSDE1, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/UpperquantileCSDE1.csv")
write.csv(LowerquantileCSDE1, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/LowerquantileCSDE1.csv")


UpperquantileCSDE1 <- read.csv(file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/UpperquantileCSDE1.csv", header = TRUE, row.names = 1)
LowerquantileCSDE1 <- read.csv(file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/LowerquantileCSDE1.csv", header = TRUE, row.names = 1)

CibersimplehighQuant <- as.data.frame(cibersimple[rownames(cibersimple) %in%  row.names(UpperquantileCSDE1),])
CibersimplelowQuant <- as.data.frame(cibersimple[rownames(cibersimple) %in%  rownames(LowerquantileCSDE1),])

write.csv(CibersimplehighQuant, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/CibersimplehighQuantcsde1high.csv")
write.csv(CibersimplelowQuant, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/CibersimplelowQuantcsde1low.csv")
