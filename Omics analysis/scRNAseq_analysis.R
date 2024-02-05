library(Seurat)
library(patchwork)
library(tidyverse)
library(ggpubr)
library(ggplot2)

#read in data
data.combined <- readRDS("/home/mjrz/Documents/Barratts_Science/20230328/data.combined.20230328.rds")


# Visualization
Idents(data.combined) <-"seurat_clusters"
DimPlot(data.combined, reduction = "umap",label = T)
DefaultAssay(data.combined) <-"RNA"
p4 <- FeaturePlot(data.combined, features =c("CSDE1"),pt.size = 1,cols =c("gold","black"),order = T, min.cutoff ="q1",max.cutoff = "q99")
p3+p4
p5 <- FeaturePlot(data.combined, features =c("HLA-A","HLA-B","HLA-C"),pt.size = 1,cols =c("gold","black"),order = T, min.cutoff ="q1",max.cutoff = "q99",split.by = "data_source")
p5

#data for HLA/CSDE1 ratio
HLAA<- FetchData(data.combined, vars = "HLA-A")
HLAB<- FetchData(data.combined, vars = "HLA-B")
HLAC<- FetchData(data.combined, vars = "HLA-C")
IRF1<- FetchData(data.combined, vars = "IRF1")
CSDE1<- FetchData(data.combined, vars = "CSDE1") 
data.combined$topHLA <- ((HLAA+0.1)+(HLAB+0.1)+(HLAC+0.1))
data.combined$dCSDE1 <- (CSDE1+0.1)
data.combined$CSDE1HLAratio <- (data.combined$dCSDE1/data.combined$topHLA)
ltest <- (data.combined$dCSDE1/data.combined$topHLA)


data.combined.3 <- subset(data.combined, subset = topHLA != "0.3")
data.combined.4 <- subset(data.combined.3, subset = dCSDE1 != "0.1")
data.combined.4$CSDE1HLAratio.norm <- (data.combined.4$CSDE1HLAratio/data.combined.4$nCount_RNA*10000)


Idents(data.combined.4) <- "merged_celltype"
v1 <- VlnPlot(data.combined.4,features = "CSDE1HLAratio.norm",group.by = "merged_celltype",idents = c( "Cancer.cells","Goblet","Columnar_Undifferentiated","Columnar_Intermediate",
                                                                                             "Columnar_differentiated","Foveolar_differentiated","Enteroendocrine_CHGA",
                                                                                             "Enteroendocrine_GHRL","Enterocytes"),log=T,sort = T)+NoLegend()+
stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95)

my_comparisons <- list(c("Foveolar_differentiated","Cancer.cells"))

v2 <-v1 + 
  stat_compare_means(comparisons = my_comparisons, na.rm=T,label = "p.signif",label.y =2 ,step.increase=0.08,vjust=0.25) +
  stat_compare_means(label.y = 1.5,label.x = "Columnar_differentiated")
ggsave(plot = v2,filename = "CSDE1_HLA_ratio_VlnPlot1.png",width = 35,height = 15,units = "cm")

data.will <- cbind(data.combined.4$CSDE1HLAratio.norm,data.combined.4$merged_celltype)
write.csv(data.will,"CSDE1HLAratio.norm.csv")

data.combined.not1  <- subset(data.combined,subset= HLAA.CSDE1.ratio != 1)
v3<-VlnPlot(data.combined.4,features = "CSDE1",group.by = "merged_celltype",idents = c( "Cancer.cells","Goblet","Columnar_Undifferentiated","Columnar_Intermediate",
                                                                                                       "Columnar_differentiated","Foveolar_differentiated","Enteroendocrine_CHGA",
                                                                                                       "Enteroendocrine_GHRL","Enterocytes"),log=T,sort = T)+NoLegend()+
  stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95)
ggsave(plot = v3,filename = "CSDE1_EAC_cells_VlnPlot1.png",width = 35,height = 15,units = "cm")

FeaturePlot(data.combined.not1, features =c("HLAA.CSDE1.ratio"),pt.size = 1,cols =c("gold","black"),order = T)

data.combined.not1$ratio.groups <- data.combined.not1$HLAA.CSDE1.ratio
data.combined.not1$ratio.groups <- cut(data.combined.not1$HLAA.CSDE1.ratio, 
                   breaks=c(-Inf,0.9,Inf), 
                   labels=c("low","high"))

DotPlot(data.combined.not1,features = c("CSDE1","HLA-A", "HLA-B","HLA-C", "IRF1","IRF9","PTPN2","INTERFERON.GAMMA.RESPONSE","INTERFERON.ALPHA.RESPONSE"),
        group.by = "ratio.groups",dot.scale = 15,idents = c( "Cancer.cells"))+theme(axis.text.x =element_text(angle=45,hjust=1, size=8))

VlnPlot(data.combined.not1,features = "IRF9",group.by = "ratio.groups",idents = c( "Cancer.cells"))+NoLegend()
VlnPlot(data.combined.not1,features = "PTPN2",group.by = "ratio.groups",idents = c( "Cancer.cells"))+NoLegend()

FeatureScatter(data.combined.not1,feature1 = "IRF1",feature2 = "HLAA.CSDE1.ratio",)
FeatureScatter(data.combined.not1,feature1 = "IRF9",feature2 = "HLAA.CSDE1.ratio")
FeatureScatter(data.combined.not1,feature1 = "PTPN2",feature2 = "HLAA.CSDE1.ratio")

DefaultAssay(data.combined) <- "RNA"
list.cells <- c(unique(metadata$IM.Status))
data.combined$IM.Status <- data.combined$IM.Status %>% replace(is.na(.), "N.A")
unique(data.combined$IM.Status)
Idents(data.combined) <-"IM.Status"

v1 <- VlnPlot(data.combined,features = "CSDE1", group.by = "IM.Status", slot="counts",sort=T,idents = c("IM with chronic gastritis","IM with HGD", "Moderately extensive IM",
                                                                                                        "No IM, Chronic gastritis & HGD", 
                                                                                                        "IM with chronic gastritis & LGD", "IM only"),log=F,y.max = 30)+NoLegend()
v1
my_comparisons <- list(c("Moderately extensive IM", "IM with chronic gastritis & LGD"), c("Moderately extensive IM", "IM with chronic gastritis"), 
                       c("Moderately extensive IM", "No IM, Chronic gastritis & HGD"))
v2<-v1 + 
  stat_compare_means(comparisons = my_comparisons, na.rm=T,label = "p.signif",label.y = 18,step.increase=0.08,vjust=0.75) +
  stat_compare_means(label.y = 25,label.x = "IM with chronic gastritis")
ggsave(plot = v2,filename = "CSDE1_VlnPlot1.png",width = 25,height = 15,units = "cm")


v1 <- VlnPlot(data.combined,features = "CSDE1", group.by = "IM.Status", slot="counts",sort=T,idents = c("IM with chronic gastritis","IM with HGD", "Moderately extensive IM",
                                                                                                        "No IM, Chronic gastritis & HGD", 
                                                                                                        "IM with chronic gastritis & LGD", "IM only"),log=F,y.max = 30)+NoLegend()
v1
my_comparisons <- list(c("IM with HGD", "IM with chronic gastritis & LGD"), c("IM with HGD", "IM with chronic gastritis"), 
                       c("IM with HGD", "No IM, Chronic gastritis & HGD"))

v2 <-v1 + 
  stat_compare_means(comparisons = my_comparisons, na.rm=T,label = "p.signif",label.y = 18,step.increase=0.08,vjust=0.75) +
  stat_compare_means(label.y = 25,label.x = "IM with chronic gastritis")
ggsave(plot = v2,filename = "CSDE1_VlnPlot2.png",width = 25,height = 15,units = "cm")

v1 <- VlnPlot(data.combined,features = "CSDE1", group.by = "IM.Status", slot="counts",sort=T,idents = c("IM with chronic gastritis","IM with HGD", "Moderately extensive IM",
                                                                                                        "No IM, Chronic gastritis & HGD", 
                                                                                                        "IM with chronic gastritis & LGD", "IM only"),log=F,y.max = 30)+NoLegend()
v1
my_comparisons <- list(c("IM only", "IM with chronic gastritis & LGD"), c("IM only", "IM with chronic gastritis"), 
                       c("IM only", "No IM, Chronic gastritis & HGD"))

v2 <-v1 + 
  stat_compare_means(comparisons = my_comparisons, na.rm=T,label = "p.signif",label.y = 18,step.increase=0.08,vjust=0.25) +
  stat_compare_means(label.y = 25,label.x = "IM with chronic gastritis")
ggsave(plot = v2,filename = "CSDE1_VlnPlot3.png",width = 25,height = 15,units = "cm")

VlnPlot(data.combined,features = "CSDE1", group.by = "IM.Status")                       
VlnPlot(data.combined.not1,features = "HLAA.CSDE1.ratio",group.by = "merged_celltype",idents = c( "Cancer","Cancer.cells","Goblet","Columnar_Undifferentiated","Columnar_Intermediate",
                                                                                             "Columnar_differentiated","Foveolar_differentiated","Enteroendocrine_CHGA",
                                                                                             "Enteroendocrine_GHRL","Enterocytes"))+NoLegend()

VlnPlot(data.combined,features = "HLAA.CSDE1.ratio",group.by = "merged_celltype",idents = c( "Cancer","Cancer.cells","Goblet","Columnar_Undifferentiated","Columnar_Intermediate",
                                                                                           "Columnar_differentiated","Foveolar_differentiated","Enteroendocrine_CHGA",
                                                                                           "Enteroendocrine_GHRL","Enterocytes"))+NoLegend()
Idents(data.combined) <-"merged_celltype"
RidgePlot(data.combined.not1, features =c("HLAA.CSDE1.ratio"),group.by = "merged_celltype",idents = c("Cancer.cells","Goblet","Columnar_Undifferentiated","Columnar_Intermediate",
                                                                                                "Columnar_differentiated","Foveolar_differentiated","Enteroendocrine_CHGA",
                                                                                                "Enteroendocrine_GHRL","Enterocytes")) +NoLegend()

VlnPlot(data.combined,features = "HLAA.CSDE1.ratio",group.by = "merged_celltype",idents = c("Cancer","Goblet","Columnar_Undifferentiated","Columnar_Intermediate",
                                                                                                    "Columnar_differentiated","Foveolar_differentiated","Enteroendocrine_CHGA",
                                                                                                    "Enteroendocrine_GHRL","Enterocytes") )+NoLegend()

VlnPlot(data.combined,features = "IFNG",group.by = "merged_celltype",idents = c("Cancer","Goblet","Columnar_Undifferentiated","Columnar_Intermediate",
                                                                                            "Columnar_differentiated","Foveolar_differentiated","Enteroendocrine_CHGA",
                                                                                            "Enteroendocrine_GHRL","Enterocytes") )+NoLegend()
VlnPlot(data.combined,features = "HMGA2",group.by = "merged_celltype",idents = c("Cancer","Goblet","Columnar_Undifferentiated","Columnar_Intermediate",
                                                                                "Columnar_differentiated","Foveolar_differentiated","Enteroendocrine_CHGA",
                                                                                "Enteroendocrine_GHRL","Enterocytes") )+NoLegend()

data.combined.flitered <- subset(data.combined, subset = merged_celltype == "Cancer.cells")
data.combined.flitered2 <- subset(data.combined.flitered, subset = data_source == "EAC")

VlnPlot(data.combined.flitered2,features = "INTERFERON.GAMMA.RESPONSE")+NoLegend()

FeatureScatter(data.combined.flitered2,feature1 = "INTERFERON.GAMMA.RESPONSE",feature2 = "HLAA.CSDE1.ratio")
FeatureScatter(data.combined.flitered2,feature1 = "INTERFERON.GAMMA.RESPONSE",feature2 = "HLA-A")

FeatureScatter(data.combined.flitered2,feature1 = "INTERFERON.ALPHA.RESPONSE",feature2 = "HLAA.CSDE1.ratio")
FeatureScatter(data.combined.flitered2,feature1 = "INTERFERON.ALPHA.RESPONSE",feature2 = "HLA-A")
FeatureScatter(data.combined.flitered2,feature1 = "INTERFERON.ALPHA.RESPONSE",feature2 = "CSDE1")

FeatureScatter(data.combined.flitered2,feature1 = "IRF1",feature2 = "CSDE1")
FeatureScatter(data.combined.flitered2,feature1 = "IRF9",feature2 = "CSDE1")
FeatureScatter(data.combined.flitered2,feature1 = "PTPN2",feature2 = "CSDE1")

data.combined.flitered <- subset(data.combined, subset = data_source == "EAC")

FeaturePlot(data.combined.flitered, features =c("HLA-A","CSDE1"),pt.size = 1,order = T, min.cutoff ="q1",max.cutoff = "q99",blend = T)
FeaturePlot(data.combined.flitered, features =c("SMYD3","CSDE1"),pt.size = 1,order = T, min.cutoff ="q1",max.cutoff = "q99",blend = T)
FeaturePlot(data.combined.flitered2, features =c("PTPN2","CSDE1"),pt.size = 1,order = T, min.cutoff ="q1",max.cutoff = "q99",blend = T)

FeaturePlot(data.combined.flitered2, features =c("HLA-A","CSDE1"),pt.size = 1,order = T, min.cutoff ="q1",max.cutoff = "q99",blend = T)


# Exploration of CSDE1 in eosphagus
DimPlot(data.combined, group.by = "merged_celltype",label = T)+NoLegend()

FeaturePlot(data.combined, features ="CSDE1", cols = c("gold", "black"),order = T,min.cutoff = "q10",max.cutoff = "q90")

FeaturePlot(data.combined, features ="CSDE1", cols = c("gold", "black"),order = T,min.cutoff = "q10",max.cutoff = "q90",split.by = "data_source")
VlnPlot(data.combined, features =c("CSDE1"),group.by = "merged_celltype",sort = T)+NoLegend()


EAC <- subset(data.combined,subset = data_source =="EAC")
VlnPlot(EAC,features = "CSDE1",group.by = "merged_celltype")

DefaultAssay(EAC) <- "RNA"
EAC <-NormalizeData(EAC)
EAC <- ScaleData(EAC)

#read in Hallmark genesets from: http://www.gsea-msigdb.org/gsea/msigdb/collections.jsp#H
hallmarks.sigs <- read.table(file="/home/mjrz/Documents/Barratts_Science/20230328/h.all.v7.5.1.symbols.txt", sep = "\t",blank.lines.skip = T,header = T)
sig.names <- colnames(hallmarks.sigs)
#Calculate per cell hallmark Module scores

for (i in sig.names) {
  module_genes <- assign(i, as.vector(hallmarks.sigs[i]))
  EAC <- AddModuleScore(object = EAC,features = module_genes,name = i, search = T)
}

for (i in sig.names) {
  colnames(EAC@meta.data) <-
    gsub(x = colnames(EAC@meta.data)
         , pattern = paste0(i,1)
         , replacement = i)
}

EAC <- AddModuleScore(object = EAC,features = INTERFERON.ALPHA.RESPONSE, name = "IFN.ALPHA.RESPONSE",search = T)
EAC <- AddModuleScore(object = EAC,features = INTERFERON.GAMMA.RESPONSE, name = "IFN.GAMMA.RESPONSE",search = T)
EAC$integrated.cell.types.manual
DefaultAssay(EAC) <- "RNA"
Idents(EAC) <- "integrated.cell.types.manual"
EAC.1 <- subset(EAC, idents = c("Cancer.cells", "T.cells"))
EAC.1 <- NormalizeData(EAC.1)
EAC.1 <- ScaleData(EAC.1)
Idents(EAC.1) <- "integrated.cell.types.manual"
EAC.2 <- AggregateExpression(EAC.1,return.seurat = T,group.by = c("orig.patient","integrated.cell.types.manual"))
EAC.2@meta.data
EAC.2$ids <- rownames(EAC.2@meta.data)

data.heatmap <- FetchData(object = EAC.2, vars = c("orig.ident", "ids","CSDE1", "PRF1","GZMB","HLA-A","HLA-B","HLA-C","IRF1","IRF9","IFNG","IFNA1",
                                                 "IFN.ALPHA.RESPONSE1","IFN.GAMMA.RESPONSE1"),layer="data")
ct <- data.frame(t(as.data.frame(strsplit(data.heatmap$ids,"_"))))
colnames(ct) <- c("orig.patient","cell.type")

data.heatmap$orig.patient <- ct$orig.patient
data.heatmap$cell.type <- ct$cell.type

write.csv(data.heatmap,file="/home/mjrz/Documents/Barratts_Science/20230328/data.heatmap.willp.csv")
Tcells.data <- subset(data.heatmap,subset=data.heatmap$integrated.cell.types.manual=="T.cells")
Cancer.data <- subset(data.heatmap,subset=data.heatmap$integrated.cell.types.manual=="Cancer.cells")


                      
annotation.1 <- data.heatmap[,1:2]
heatmap.matrix <- as.matrix(t(Tcells.data[,4:13]))
library(pheatmap)
#library(RColorBrewer)

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

data_subset_norm <- t(apply(heatmap.matrix, 1, cal_z_score))

# Define color scale using quantiles
my_breaks <- quantile(as.matrix(data_subset_norm),
                      probs = seq(0, 1, by = 0.1))
my_breaks #trim
my_breaks2 <- c(-2, -0.8347214 ,-0.5512894 ,-0.3581906, -0.2340996 ,-0.2243145 ,-0.1336555,  0.1948927,  0.6218952,  1.0867393, 2)
my_palette <- colorRampPalette(c("blue","white","red"))(11)

pheatmap(data_subset_norm,show_colnames = FALSE,annotation_col = annotation.1,scale = "none",color = my_palette,breaks = my_breaks2,cluster_cols = FALSE)
pheatmap(data_subset_norm,show_colnames = FALSE,annotation_col = annotation.1,scale = "none",color = my_palette,
         breaks = my_breaks2,cluster_cols = F,cluster_rows = F)

#####
annotation.1 <- data.heatmap[,1:2]
heatmap.matrix <- as.matrix(t(Cancer.data[,c(3,6,7,8)]))
library(pheatmap)
#library(RColorBrewer)

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

data_subset_norm <- t(apply(heatmap.matrix, 1, cal_z_score))

# Define color scale using quantiles
my_breaks <- quantile(as.matrix(data_subset_norm),
                      probs = seq(0, 1, by = 0.1))
my_breaks #trim
my_breaks2 <- c(-2, -0.8347214 ,-0.5512894 ,-0.3581906, -0.2340996 ,-0.2243145 ,-0.1336555,  0.1948927,  0.6218952,  1.0867393, 2)
my_palette <- colorRampPalette(c("blue","white","red"))(11)

pheatmap(data_subset_norm,show_colnames = FALSE,annotation_col = annotation.1,scale = "none",color = my_palette,breaks = my_breaks2,cluster_cols = F)

#####
data.heatmap.sort <- data.heatmap[order(data.heatmap$CSDE1),]
Tcells.data <- subset(data.heatmap.sort,subset=data.heatmap.sort$cell.type=="T.cells")
Cancer.data$orig.ident <- factor(Cancer.data$orig.ident,levels = c(Cancer.data$orig.ident))
levels(Cancer.data$orig.ident)
c(Cancer.data$orig.ident)
levels(Tcells.data$orig.ident)
c(Tcells.data$orig.ident)
Tcells.data <- Tcells.data[Tcells.data$orig.ident %in% Cancer.data$orig.ident,]
Tcells.data$orig.ident <- factor(Tcells.data$orig.ident,levels = c(Tcells.data$orig.ident))
Tcells.data <- Tcells.data[match(levels(Cancer.data$orig.ident),Tcells.data$orig.ident),]
Tcells.data$orig.ident <- factor(Tcells.data$orig.ident,levels = c(Tcells.data$orig.ident))
annotation.1 <- Tcells.data[,c(1,13)]
heatmap.matrix <- as.matrix(t(Tcells.data[,c(4,5,9,10,11)]))

library(pheatmap)
#library(RColorBrewer)

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

data_subset_norm <- t(apply(heatmap.matrix, 1, cal_z_score))

# Define color scale using quantiles
my_breaks <- quantile(as.matrix(data_subset_norm),
                      probs = seq(0, 1, by = 0.1))
my_breaks #trim
my_breaks2 <- c(-2, -1.0348679, -0.9271195, -0.6786898 ,-0.1330739,  0.1932113 , 0.4443534,  0.6004274,  0.7706185,  1.0137190,  2)
my_palette <- colorRampPalette(c("blue","white","red"))(11)

pheatmap(data_subset_norm,show_colnames = T,annotation_col = annotation.1,scale = "none",color = my_palette,breaks = my_breaks2,cluster_rows = T,
         cluster_cols = F)

annotation.1 <- Cancer.data[,c(1,13)]
heatmap.matrix <- as.matrix(t(Cancer.data[,c(3,6,7,8)]))
library(pheatmap)
#library(RColorBrewer)

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

data_subset_norm <- t(apply(heatmap.matrix, 1, cal_z_score))

# Define color scale using quantiles
my_breaks <- quantile(as.matrix(data_subset_norm),
                      probs = seq(0, 1, by = 0.1))
my_breaks #trim
my_breaks2 <- c(-2, -1.0348679, -0.9271195, -0.6786898 ,-0.1330739,  0.1932113 , 0.4443534,  0.6004274,  0.7706185,  1.0137190,  2)
my_palette <- colorRampPalette(c("blue","white","red"))(11)

pheatmap(data_subset_norm,show_colnames = T,annotation_col = annotation.1,scale = "none",color = my_palette,breaks = my_breaks2,
         cluster_rows = F,
         cluster_cols = F)
