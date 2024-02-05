allmarkersHscore <- read.csv("C:/Users/wp1g19/Documents/Revisit bioinformatics/TMA combined heatmap/ALLMARKERSHscore.csv", header = TRUE, row.names = 1)
allmarkersHscore <- na.omit(allmarkersHscore)
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

allmarkersHscoreZ <- t(apply(allmarkersHscore, 1, cal_z_score))
allmarkersHscoreZ <- t(allmarkersHscoreZ)
my_breaks <- quantile(as.matrix(allmarkersHscoreZ),
                      probs = seq(0, 1, by = 0.1))
my_breaks #trim
my_breaks2 <- c(-3, -0.76171131, -0.65291548, -0.57996258, -0.52375073, -0.44444683, -0.28360139, -0.01150961,  0.72719666,  1.71735247,  3)
my_palette <- colorRampPalette(c("blue","white","red"))(10)

out <- pheatmap(allmarkersHscoreZ,treeheight_col = 50, fontsize = 15,breaks = my_breaks2, show_colnames = F,cluster_rows = T, cluster_cols = T,scale = "row",color = my_palette, clustering_method = "ward.D2",clustering_distance_rows = "canberra",clustering_distance_cols = "canberra" ,main = "OAC cell line MHC class I/regulator gene expression profile")
dev.off()

clusters <- cutree(out$tree_col, k=4)[out$tree_col[["order"]]]
annot_col <- data.frame(col.names = names(clusters),
                        cluster = as.factor(clusters))
annot_col$cluster <- gsub("3","5", annot_col$cluster)
annot_col$cluster <- gsub("1","6", annot_col$cluster)
annot_col$cluster <- gsub("4","7", annot_col$cluster)
annot_col$cluster <- gsub("2","8", annot_col$cluster)
annot_col$cluster <- gsub("5","1", annot_col$cluster)
annot_col$cluster <- gsub("6","2", annot_col$cluster)
annot_col$cluster <- gsub("7","3", annot_col$cluster)
annot_col$cluster <- gsub("8","4", annot_col$cluster)
annot_col$col.names <- NULL

annot_colors <- list(cluster=c("1"="blue","2"="green","3"="red","4"="purple"))

pheatmap(allmarkersHscoreZ,treeheight_col = 50,annotation_col =annot_col,annotation_colors = annot_colors,  fontsize = 15,breaks = my_breaks2, show_colnames = F,cluster_rows = T, cluster_cols = T,scale = "row",color = my_palette, clustering_method = "ward.D2",clustering_distance_rows = "canberra",clustering_distance_cols = "canberra" ,main = "APM H scores and T cell composition (n=78)")
