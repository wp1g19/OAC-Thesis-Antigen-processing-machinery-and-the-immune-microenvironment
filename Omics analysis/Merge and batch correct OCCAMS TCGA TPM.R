# Gather TCGA data
library(TCGAbiolinks)
library(SummarizedExperiment)
#Filter to EAC
setwd("C:/Users/wp1g19/Documents/Revisit bioinformatics/data")
esca.subtype_OAC <- TCGAquery_subtype("esca")
subset_OAC_filter <- subset(esca.subtype_OAC, esca.subtype_OAC$`Histological Type - Oesophagus` == 'EAC')
query.exp <- GDCquery(project = "TCGA-ESCA", 
                      legacy = FALSE,
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "STAR - Counts", sample.type = "Primary Tumor",
                      barcode =  c("TCGA-L7-A6VZ", "TCGA-IG-A4QS", "TCGA-L5-A8NV", "TCGA-R6-A6XG", "TCGA-R6-A6DQ", "TCGA-L5-A8NE", "TCGA-L5-A8NT", "TCGA-L5-A4OQ", "TCGA-2H-A9GH", "TCGA-R6-A8W8", "TCGA-L5-A8NJ", "TCGA-V5-AASX", "TCGA-L5-A4OP", "TCGA-L5-A8NF", "TCGA-L5-A8NM", "TCGA-L5-A8NR", "TCGA-2H-A9GN", "TCGA-L5-A8NS", "TCGA-L5-A4OI", "TCGA-2H-A9GO", "TCGA-R6-A8W5", "TCGA-L5-A88V", "TCGA-L5-A8NI", "TCGA-L5-A4OW", "TCGA-L5-A4OJ", "TCGA-L5-A8NW", "TCGA-R6-A6L6", "TCGA-L5-A891", "TCGA-RE-A7BO", "TCGA-IC-A6RE", "TCGA-JY-A6F8", "TCGA-JY-A93E", "TCGA-VR-AA4D", "TCGA-R6-A6KZ", "TCGA-S8-A6BV", "TCGA-X8-AAAR", "TCGA-L5-A4OO", "TCGA-JY-A6FB", "TCGA-JY-A93C", "TCGA-2H-A9GK", "TCGA-2H-A9GF", "TCGA-2H-A9GI", "TCGA-2H-A9GL", "TCGA-2H-A9GM", "TCGA-2H-A9GR","TCGA-JY-A6FH", "TCGA-JY-A93D", "TCGA-JY-A938", "TCGA-L5-A4OE", "TCGA-L5-A4OF", "TCGA-L5-A4OG","TCGA-L5-A4OH", "TCGA-L5-A4ON", "TCGA-L5-A4OR", "TCGA-L5-A4OS", "TCGA-L5-A4OT", "TCGA-L5-A4OU", "TCGA-L5-A4OX", "TCGA-L5-A8NG", "TCGA-L5-A8NH", "TCGA-L5-A8NL", "TCGA-L5-A8NN", "TCGA-L5-A8NU", "TCGA-L5-A43C", "TCGA-L5-A43I", "TCGA-L5-A88Y", "TCGA-M9-A5M8", "TCGA-Q9-A6FW", "TCGA-R6-A6DN", "TCGA-R6-A6L4", "TCGA-R6-A6XQ", "TCGA-R6-A6Y2", "TCGA-R6-A8WC", "TCGA-R6-A8WG", "TCGA-V5-A7RB", "TCGA-V5-A7RE", "TCGA-V5-AASW", "TCGA-VR-A8EQ", "TCGA-ZR-A9CJ"))
GDCdownload(query.exp)
OACESCA.exp <- GDCprepare(query = query.exp, save = TRUE, save.filename = "OACESCAExp.rda")
library(dplyr)
OACESCA.exp2 <- assay(OACESCA.exp,"unstranded") 

dataPrep <- TCGAanalyze_Preprocessing(object = OACESCA.exp)
#convert gene IDs
library(biomaRt)
dataPrep <- as.data.frame(dataPrep)
dataPrep$ensembl <- gsub("\\..*","", rownames(dataPrep))
dataPrep = dataPrep[!duplicated(dataPrep$ensembl),]
rownames(dataPrep) <- dataPrep$ensembl
dataPrep$ensembl <- NULL
mart <- useEnsembl(biomart = "ensembl", 
                   dataset = "hsapiens_gene_ensembl", 
                   mirror = "useast")
genes <- rownames(dataPrep)
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
G_list <- G_list[!(is.na(G_list$hgnc_symbol) | G_list$hgnc_symbol==""), ]
G_list <- G_list[!duplicated(G_list$ensembl_gene_id),]
combined_df_removedcols2 <- dataPrep[ rownames(dataPrep) %in% G_list$ensembl_gene_id ,]
combined_df_removedcols3 <- subset(combined_df_removedcols2, rownames(combined_df_removedcols2) != "ENSG00000230417")
genes <- rownames(combined_df_removedcols3)

G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
G_list[ duplicated(G_list$ensembl_gene_id) | duplicated(G_list$ensembl_gene_id, fromLast = TRUE), ]
G_list <- distinct(G_list, ensembl_gene_id, .keep_all = TRUE) 
rownames(G_list) <- G_list$ensembl_gene_id
combined_df_removedcols <- merge(G_list, combined_df_removedcols3, by = 0)
#remove non-unique gene IDs
combined_df_removedcols = combined_df_removedcols[!duplicated(combined_df_removedcols$hgnc_symbol),]
rownames(combined_df_removedcols) <- combined_df_removedcols$hgnc_symbol
library(edgeR)
# make data frame and numeric
library(dplyr)
combined_df_removedcols$Row.names <- NULL
combined_df_removedcols$ensembl_gene_id <- NULL
combined_df_removedcols$hgnc_symbol <- NULL
combined_df_removedcols$NA. <- NULL
combined_df_removedcols$NA..1 <- NULL
combined_df_removedcols$NA..2 <- NULL
combined_df_removedcols$NA..3 <- NULL
combined_df_removedcols <- mutate_all(combined_df_removedcols, function(x) as.numeric(as.character(x)))

library(edgeR)
dge <- DGEList(counts = combined_df_removedcols)
keep  <- filterByExpr(dge)
table(keep)
dge <- dge[keep, , keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge,method = "TPM", lib.size=NULL)
TCGATMMMatrix <- cpm(dge)


# OCCAMS data
mypath = "C:/Users/wp1g19/Desktop/OCAAMS_RNA"

setwd(mypath)

# Create list of text files
txt_files_ls = list.files(path=mypath, pattern="*.txt") 

# Read the files in, assuming tab delimited
txt_files_df <- lapply(txt_files_ls, function(x) {read.table(file = x, header = T, sep ="\t")})
# Combine them
combined_df <- do.call("cbind", lapply(txt_files_df, as.data.frame))

combined_list_removedcols <- lapply(txt_files_df, function(x) { x[c("readcounts_intersectionNotEmpt", "genelength", "R.FPKM_union", "R.FPKM_intersectionNotEmpt")] <- NULL; x })

combined_df_removedcols <- do.call("cbind", lapply(combined_list_removedcols, as.data.frame))

colnames(combined_df_removedcols) <- txt_files_ls

#convert gene IDs
library(biomaRt)
combined_df_removedcols <- as.data.frame(combined_df_removedcols)
rownames(combined_df_removedcols) <- gsub("\\..*","", rownames(combined_df_removedcols))
combined_df_removedcols = combined_df_removedcols[!duplicated(rownames(combined_df_removedcols)),]
mart <- useEnsembl(biomart = "ensembl", 
                   dataset = "hsapiens_gene_ensembl", 
                   mirror = "useast")
genes <- rownames(combined_df_removedcols)
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
G_list <- G_list[!(is.na(G_list$hgnc_symbol) | G_list$hgnc_symbol==""), ]
G_list <- G_list[!duplicated(G_list$ensembl_gene_id),]
combined_df_removedcols2 <- combined_df_removedcols[ rownames(combined_df_removedcols) %in% G_list$ensembl_gene_id ,]
combined_df_removedcols3 <- subset(combined_df_removedcols2, rownames(combined_df_removedcols2) != "ENSG00000230417")
genes <- rownames(combined_df_removedcols3)
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
rownames(G_list) <- G_list$ensembl_gene_id
combined_df_removedcols <- merge(G_list, combined_df_removedcols3, by = 0)
#remove non-unique gene IDs
combined_df_removedcols = combined_df_removedcols[!duplicated(combined_df_removedcols$hgnc_symbol),]
rownames(combined_df_removedcols) <- combined_df_removedcols$hgnc_symbol
library(edgeR)
# make data frame and numeric
library(dplyr)
combined_df_removedcols$Row.names <- NULL
combined_df_removedcols$ensembl_gene_id <- NULL
combined_df_removedcols$hgnc_symbol <- NULL
combined_df_removedcols$NA. <- NULL
combined_df_removedcols$NA..1 <- NULL
combined_df_removedcols$NA..2 <- NULL
combined_df_removedcols$NA..3 <- NULL
combined_df_removedcols <- mutate_all(combined_df_removedcols, function(x) as.numeric(as.character(x)))

OCCAMSMatrixnoNA <- combined_df_removedcols

#/ make the DGEList:

y <- DGEList(OCCAMSMatrixnoNA)

keep <- filterByExpr(y)

y <- y[keep, , keep.lib.sizes=FALSE]
#/ calculate TMM normalization factors:
y <- calcNormFactors(y, method = "TPM")

#/ get the normalized counts:
OCCAMSTPMMatrixnoNA <- cpm(y, log=FALSE)

colnames(OCCAMSTPMMatrixnoNA) <- gsub("_1.0_rpkm.txt*","", colnames(OCCAMSTPMMatrixnoNA))

# Convert the SLX to OCCAMS number
OCCAMSTPMMatrixnoNA2 <- t(OCCAMSTPMMatrixnoNA)
conversionmap <- read.csv("C:/Users/wp1g19/Desktop/Post return writing/TPR/Cibersortx/SLXtoIDOCCAMS.csv")
rownames(conversionmap) <- conversionmap$SLX_ID
rownames(conversionmap) <- gsub("-",".", rownames(conversionmap))
rownames(OCCAMSTPMMatrixnoNA2) <- gsub("-",".", rownames(OCCAMSTPMMatrixnoNA2))
rownames(conversionmap) <- gsub("_",".", rownames(conversionmap))
rownames(OCCAMSTPMMatrixnoNA2) <- gsub("_",".", rownames(OCCAMSTPMMatrixnoNA2))

OCCAMSTPMMatrixnoNA3 <- as.data.frame(OCCAMSTPMMatrixnoNA2)
OCCAMSTPMMatrixnoNA3 <- OCCAMSTPMMatrixnoNA2[ rownames(OCCAMSTPMMatrixnoNA2) %in% rownames(conversionmap),]
conversionmaporder <- conversionmap[rownames(conversionmap) %in% row.names(OCCAMSTPMMatrixnoNA3) ,]

row.names(OCCAMSTPMMatrixnoNA3) <- conversionmaporder$ID
OCCAMSTMMMatrixnoNA <- as.data.frame(t(OCCAMSTPMMatrixnoNA3))

batch <- c(rep("A",79), rep("B", 141))

library(pcaExplorer)
library(mixOmics)
library(pcaExplorer)
library(PLSDAbatch)
pca <- mixOmics::pca(t(OCCAMSTMMMatrixnoNA), ncomp = 3)


Scatter_Density(pca)#

outlierspcacheck <- as.data.frame(pca$variates)
#Highest PC1 sample OC/NT/193 OC/AH/200 OC/RS/032 OC/AH/197 OC/QE/148
#Lowest 5 PC2 samples OC/AH/234 OC/AH/204 OC/AH/279 OC/RS/096 OC/SH/170

outlier_barcode <-c("OC/AH/408",
                    "OC/NT/193",
                    "OC/AH/200",
                    "OC/RS/032",
                    "OC/AH/197",
                    "OC/AH/197",
                    "OC/NT/190",
                    "OC/AH/234",
                    "OC/AH/204",
                    "OC/AH/279",
                    "OC/RS/094",
                    "OC/RS/031",
                    "OC/QE/111",
                    "OC/SH/159",
                    "OC/AH/255",
                    "OC/QE/173",
                    "OC/ED/163",
                    "OC/AH/127",
                    "OC/AH/180",
                    "OC/AH/218",
                    "OC/AH/280",
                    "OC/ED/010",
                    "OC/ED/105",
                    "OC/ED/139",
                    "OC/NS/019",
                    "OC/NT/126",
                    "OC/RS/024",
                    "OC/RS/036",
                    "OC/RS/050",
                    "OC/RS/076",
                    "OC/SH/051",
                    "OC/SH/058",
                    "OC/SH/137",
                    "OC/SH/175",
                    "OC/ST/016",
                    "OC/ST/044",
                    "OC/ST/070"
                    )

OCCAMSTMMMatrixnoNAoutliercut <- OCCAMSTMMMatrixnoNA[,!colnames(OCCAMSTMMMatrixnoNA) %in% outlier_barcode ]

OCCAMSsites <- as.data.frame(colnames(OCCAMSTMMMatrixnoNAoutliercut))
colnames(OCCAMSsites) <- "OCCAMSID"
OCCAMSsites$sitecode <- as.factor(strtrim(OCCAMSsites$OCCAMSID, 5))
summary(OCCAMSsites$sitecode)
batch <- c(rep("AH",41),
           rep("CO", 1),
           rep("ED", 5),
           rep("GS", 3),
           rep("NN", 1),
           rep("NT", 10),
           rep("PL", 1),
           rep("QE", 3),
           rep("RS", 20),
           rep("SA", 4),
           rep("SH", 10),
           rep("ST", 6),
           rep("UC", 5),
           rep("WG", 3),
           rep("WO", 1),
           rep("WY", 1))


pca <- mixOmics::pca(t(OCCAMSTMMMatrixnoNAoutliercut), ncomp = 3)
colorlist <- rainbow(17)
Scatter_Density(pca, batch = batch, color.set = colorlist)
#After removing poor quality samples no batch effect between OCCAMS sites was discovered.

#correlate expression of CSDE1 and APM genes in the occams data set

cor1 <- as.matrix(t(OCCAMSTMMMatrixnoNAoutliercut))

csdgenes <- c("CSDE1", "RXF5", "RFXAP", "RFXANK", "NLRC5", "IRF1", "CIITA", "HLA-A", "HLA-B", "HLA-C", "HLA-E", "TAP1", "TAP2", "TAPBP", "B2M")

cor1 <- OCCAMSTMMMatrixnoNAoutliercut[rownames(OCCAMSTMMMatrixnoNAoutliercut) %in% csdgenes, ]
cor_mat <- cor(t(cor1))                # Correlation matrix of example data
cor_mat                              # Print correlation matrix


library("psych")                     # Load psych package

cor_test_mat <- corr.test(t(cor1))$p    # Apply corr.test function
cor_test_mat                         # Print matrix of p-values


library("corrplot")                  # Load corrplot package
corrplot(cor_mat)                    # Draw corrplot
corrplot(cor_mat,                    # Draw corrplot with p-values
         p.mat = cor_test_mat,
         insig = "p-value")


library("ggcorrplot")                # Load ggcorrplot package

ggcorrplot(cor_mat)                  # Draw ggcorrplot
while (dev.cur()>1) dev.off()
pdf("Combinedexpressioncorrellationheatmap.pdf", height = 11, width = 8.5)
p1 <- ggcorrplot(cor_mat,                  # Draw ggcorrplot with p-values
                 p.mat = cor_test_mat, show.legend = TRUE, sig.level = 0.05, lab = TRUE,insig = "pch", pch = 4,pch.cex = 2, lab_size = 0.5)
while (dev.cur()>1) dev.off()
p2 <- ggcorrplot(cor_mat,                  # Draw ggcorrplot with p-values
                 p.mat = cor_test_mat, type ="lower",show.legend = TRUE, lab = TRUE, sig.level = 0.05, pch = 0,pch.cex = 0, lab_size = 5)

library(ggpubr)
regOC <- as.data.frame(t(OCCAMSTMMMatrixnoNAoutliercut))

quantile(regOC$`HLA-A`, probs = seq(.1, .9, by = .1))
quantile(regOC$CSDE1, probs = seq(.1, .9, by = .1))
regOC1 <- regOC[regOC$CSDE1 >50,]
regOC2 <- regOC1[regOC1$CSDE1 <1000,]
regOC3 <- regOC2[regOC2$`HLA-A` <3000,]
regOC4 <- regOC3[regOC3$`HLA-A` >300,]

ggplot(regOC4, aes(x=CSDE1, y=`HLA-A`)) + 
  geom_point() +stat_smooth(method = "lm", col = "red") +
  stat_cor(label.x = 3, label.y = 2000) 

library(mixOmics)
#correlate expression of CSDE1 and APM genes in the TCGA data set
TCGATMMMatrix <- as.data.frame(TCGATMMMatrix)
pca <- mixOmics::pca(TCGATMMMatrix, ncomp = 3)
colorlist <- rainbow(17)
Scatter_Density(pca)


cor1 <- as.matrix(t(TCGATMMMatrix))

#csdgenes <- c("HLA-A", "HLA-B", "HLA-C", "HLA-E", "TAP1", "TAP2", "TAPBP", "B2M", "CSDE1")

cor1 <- TCGATMMMatrix[rownames(TCGATMMMatrix) %in% csdgenes, ]
cor_mat <- cor(t(cor1))                # Correlation matrix of example data
cor_mat                              # Print correlation matrix


library("psych")                     # Load psych package

cor_test_mat <- corr.test(t(cor1))$p    # Apply corr.test function
cor_test_mat                         # Print matrix of p-values


library("corrplot")                  # Load corrplot package
corrplot(cor_mat)                    # Draw corrplot
corrplot(cor_mat,                    # Draw corrplot with p-values
         p.mat = cor_test_mat,
         insig = "p-value")


library("ggcorrplot")                # Load ggcorrplot package

ggcorrplot(cor_mat)                  # Draw ggcorrplot
while (dev.cur()>1) dev.off()
pdf("TCGAexpressioncorrellationheatmap.pdf", height = 11, width = 8.5)
p1 <- ggcorrplot(cor_mat,                  # Draw ggcorrplot with p-values
                 p.mat = cor_test_mat, show.legend = TRUE, sig.level = 0.05, lab = TRUE,insig = "pch", pch = 4,pch.cex = 2, lab_size = 0.5)
while (dev.cur()>1) dev.off()
p2 <- ggcorrplot(cor_mat,                  # Draw ggcorrplot with p-values
                 p.mat = cor_test_mat, type ="lower",show.legend = TRUE, lab = TRUE, sig.level = 0.05, pch = 0,pch.cex = 0, lab_size = 5)

library(ggpubr)
regOC <- as.data.frame(t(TCGATMMMatrix))

quantile(regOC$`HLA-A`, probs = seq(.1, .9, by = .1))
quantile(regOC$CSDE1, probs = seq(.1, .9, by = .1))
regOC1 <- regOC[regOC$CSDE1 >50,]
regOC2 <- regOC1[regOC1$CSDE1 <1000,]
regOC3 <- regOC2[regOC2$`HLA-A` <1500,]
regOC4 <- regOC3[regOC3$`HLA-A` >300,]

ggplot(regOC2, aes(x=CSDE1, y=`HLA-A`)) + 
  geom_point() +stat_smooth(method = "lm", col = "red") +
  stat_cor(label.x = 3, label.y = 2000) 

#join OCCAMS and TCGA and batch correct
mergeTCGA <- as.data.frame(TCGATMMMatrix)
mergeOCCAMS <- as.data.frame(OCCAMSTMMMatrixnoNAoutliercut)

mergeTCGA$GeneSymbol <- rownames(mergeTCGA)
mergeOCCAMS$GeneSymbol <- rownames(mergeOCCAMS)


TCGAOCCAMSmerged <- left_join(mergeTCGA, mergeOCCAMS, by = "GeneSymbol")
rownames(TCGAOCCAMSmerged) <- TCGAOCCAMSmerged$GeneSymbol
TCGAOCCAMSmerged$GeneSymbol <- NULL
TCGAOCCAMSmerged <- na.omit(TCGAOCCAMSmerged)

Datasets <- as.data.frame(colnames(TCGAOCCAMSmerged))
colnames(Datasets) <- "Sample ID"
Datasets$Dataset <- as.factor(strtrim(Datasets$`Sample ID`, 5))

batch <- c(rep("TCGA",79),
           rep("OCCAMS", 115))


pca <- mixOmics::pca(t(TCGAOCCAMSmerged), ncomp = 3)
colorlist <- rainbow(2)
Scatter_Density(pca, batch = batch, color.set = colorlist)

# Batch correct the 2 datasets and pca

TCGAOCCAMSmergedLimmaBatchC <- limma::removeBatchEffect(TCGAOCCAMSmerged, batch = batch)
pca <- mixOmics::pca(t(TCGAOCCAMSmergedLimmaBatchC), ncomp = 3)
colorlist <- rainbow(2)
Scatter_Density(pca, batch = batch, color.set = colorlist)

library(sva)

TCGAOCCAMSmergedComBatBatchC <- sva::ComBat(TCGAOCCAMSmerged,
                                            batch = batch,
                                            prior.plots = TRUE,
                                            mean.only = FALSE)

pca <- mixOmics::pca(t(TCGAOCCAMSmergedComBatBatchC), ncomp = 3)
colorlist <- rainbow(2)
Scatter_Density(pca, batch = batch, color.set = colorlist)

write.csv(TCGAOCCAMSmergedComBatBatchC, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/OCCAMSTCGATPMMergedBatchCorrected.csv")
#assess the CSDE1 APM correlation in the combined batch corrected data

cor1 <- as.matrix(TCGAOCCAMSmergedComBatBatchC)

csdgenes <- c("CSDE1", "RXF5", "RFXAP", "RFXANK", "NLRC5", "IRF1", "CIITA", "HLA-A", "HLA-B", "HLA-C", "HLA-E", "TAP1", "TAP2", "TAPBP", "B2M")

cor1 <- as.data.frame(TCGAOCCAMSmergedComBatBatchC[rownames(TCGAOCCAMSmergedComBatBatchC) %in% csdgenes, ])
cor_mat <- cor(t(cor1))                # Correlation matrix of example data
cor_mat                              # Print correlation matrix


library("psych")                     # Load psych package

cor_test_mat <- corr.test(t(cor1))$p    # Apply corr.test function
cor_test_mat                         # Print matrix of p-values


library("corrplot")                  # Load corrplot package
corrplot(cor_mat)                    # Draw corrplot
corrplot(cor_mat,                    # Draw corrplot with p-values
         p.mat = cor_test_mat,
         insig = "p-value")


library("ggcorrplot")                # Load ggcorrplot package

ggcorrplot(cor_mat)                  # Draw ggcorrplot
while (dev.cur()>1) dev.off()
pdf("combinedexpressioncorrellationheatmap.pdf", height = 11, width = 8.5)
p1 <- ggcorrplot(cor_mat,                  # Draw ggcorrplot with p-values
                 p.mat = cor_test_mat, show.legend = TRUE, sig.level = 0.05, lab = TRUE,insig = "pch", pch = 4,pch.cex = 2, lab_size = 0.5)
while (dev.cur()>1) dev.off()
p2 <- ggcorrplot(cor_mat,                  # Draw ggcorrplot with p-values
                 p.mat = cor_test_mat, type ="lower",show.legend = TRUE, lab = TRUE, sig.level = 0.05, pch = 0,pch.cex = 0, lab_size = 5)

library(ggpubr)
regOC <- as.data.frame(t(TCGAOCCAMSmergedComBatBatchC))

quantile(regOC$`HLA-A`, probs = seq(.1, .9, by = .1))
quantile(regOC$CSDE1, probs = seq(.1, .9, by = .1))
regOC1 <- regOC[regOC$CSDE1 >50,]
regOC2 <- regOC1[regOC1$CSDE1 <1000,]
regOC3 <- regOC2[regOC2$`HLA-A` <2000,]
regOC4 <- regOC3[regOC3$`HLA-A` >1,]

ggplot(regOC4, aes(x=CSDE1, y=`HLA-A`)) + 
  geom_point() +stat_smooth(method = "lm", col = "red") +
  stat_cor(label.x = 3, label.y = 1500) 

cor1 <- as.data.frame(t(regOC4))
cor1 <- as.data.frame(cor1[rownames(cor1) %in% csdgenes, ])
cor_mat <- cor(t(cor1))                # Correlation matrix of example data
cor_mat                              # Print correlation matrix


library("psych")                     # Load psych package

cor_test_mat <- corr.test(t(cor1))$p    # Apply corr.test function
cor_test_mat                         # Print matrix of p-values


library("corrplot")                  # Load corrplot package
corrplot(cor_mat)                    # Draw corrplot
corrplot(cor_mat,                    # Draw corrplot with p-values
         p.mat = cor_test_mat,
         insig = "p-value")


library("ggcorrplot")                # Load ggcorrplot package

ggcorrplot(cor_mat)                  # Draw ggcorrplot
while (dev.cur()>1) dev.off()
pdf("combinedexpressioncorrellationheatmapoutlierscut.pdf", height = 11, width = 8.5)
p1 <- ggcorrplot(cor_mat,                  # Draw ggcorrplot with p-values
                 p.mat = cor_test_mat, show.legend = TRUE, sig.level = 0.05, lab = TRUE,insig = "pch", pch = 4,pch.cex = 2, lab_size = 0.5)
while (dev.cur()>1) dev.off()
p2 <- ggcorrplot(cor_mat,                  # Draw ggcorrplot with p-values
                 p.mat = cor_test_mat, type ="lower",show.legend = TRUE, lab = TRUE, sig.level = 0.05, pch = 0,pch.cex = 0, lab_size = 5)

regcut <- as.data.frame(t(cor1))
regcut1 <- regcut[regcut$CSDE1 >50,]
regcut2 <- regcut1[regcut1$CSDE1 <1000,]
regcut3 <- regcut2[regcut2$`HLA-A` <2000,]
regcut4 <- regcut3[regcut3$`HLA-A` >1,]
regcut5 <- regcut4[regcut4$`HLA-B` <4000,]
regcut6 <- regcut5[regcut5$`HLA-C` <4000,]

ggplot(regcut6, aes(x=CSDE1, y=)) + 
  geom_point() +stat_smooth(method = "lm", col = "red") +
  stat_cor(label.x = 3, label.y = 1500) 

cor1 <- as.matrix(t(regcut6))


cor_mat <- cor(t(cor1))                # Correlation matrix of example data
cor_mat                           # Print correlation matrix

#order cor mat

library("psych")                     # Load psych package

cor_test_mat <- corr.test(t(cor1))$p    # Apply corr.test function
cor_test_mat             # Print matrix of p-values



library("corrplot")                  # Load corrplot package
corrplot(cor_mat)                    # Draw corrplot
corrplot(cor_mat,                    # Draw corrplot with p-values
         p.mat = cor_test_mat,
         insig = "p-value")


library("ggcorrplot")                # Load ggcorrplot package

ggcorrplot(cor_mat)                  # Draw ggcorrplot
while (dev.cur()>1) dev.off()
pdf("combinedexpressioncorrellationheatmapOutliersCut.pdf", height = 11, width = 8.5)
p1 <- ggcorrplot(cor_mat,                  # Draw ggcorrplot with p-values
                 p.mat = cor_test_mat, show.legend = TRUE, sig.level = 0.05, lab = FALSE,insig = "pch", pch = 4,pch.cex = 3, lab_size = 2)
while (dev.cur()>1) dev.off()
p2 <- ggcorrplot(cor_mat,                  # Draw ggcorrplot with p-values
                 p.mat = cor_test_mat, type ="upper",show.legend = TRUE, lab = TRUE, sig.level = 0.05, pch = 1,pch.cex = 0,pch.col = "red", lab_size = 5, hc.order = FALSE)
p2 + ggtitle("Correlation of MHC I genes and CSDE1 n = 176")
regcut6$NLRC5
ggplot(regcut6, aes(x=CSDE1, y=IRF1)) + 
  geom_point() +stat_smooth(method = "lm", col = "red") +
  stat_cor(label.x = 3, label.y = 1500) 

#quick rewrite of ggcorplot to label the significant values not the insignificant values

# function body
ggcorrplot <- function(corr,
                       method = c("square", "circle"),
                       type = c("full", "lower", "upper"),
                       ggtheme = ggplot2::theme_minimal,
                       title = "",
                       show.legend = TRUE,
                       legend.title = "Corr",
                       show.diag = NULL,
                       colors = c("blue", "white", "red"),
                       outline.color = "gray",
                       hc.order = FALSE,
                       hc.method = "complete",
                       lab = FALSE,
                       lab_col = "black",
                       lab_size = 4,
                       p.mat = NULL,
                       sig.level = 0.05,
                       insig = c("pch", "blank"),
                       pch = 4,
                       pch.col = "black",
                       pch.cex = 5,
                       tl.cex = 12,
                       tl.col = "black",
                       tl.srt = 45,
                       digits = 2,
                       as.is = FALSE) {
  type <- match.arg(type)
  method <- match.arg(method)
  insig <- match.arg(insig)
  if (is.null(show.diag)) {
    if (type == "full") {
      show.diag <- TRUE
    } else {
      show.diag <- FALSE
    }
  }
  
  if (inherits(corr, "cor_mat")) {
    # cor_mat object from rstatix
    cor.mat <- corr
    corr <- .tibble_to_matrix(cor.mat)
    p.mat <- .tibble_to_matrix(attr(cor.mat, "pvalue"))
  }
  
  if (!is.matrix(corr) & !is.data.frame(corr)) {
    stop("Need a matrix or data frame!")
  }
  corr <- as.matrix(corr)
  
  corr <- base::round(x = corr, digits = digits)
  
  if (hc.order) {
    ord <- .hc_cormat_order(corr, hc.method = hc.method)
    corr <- corr[ord, ord]
    if (!is.null(p.mat)) {
      p.mat <- p.mat[ord, ord]
      p.mat <- base::round(x = p.mat, digits = digits)
    }
  }
  
  if (!show.diag) {
    corr <- .remove_diag(corr)
    p.mat <- .remove_diag(p.mat)
  }
  
  # Get lower or upper triangle
  if (type == "lower") {
    corr <- .get_lower_tri(corr, show.diag)
    p.mat <- .get_lower_tri(p.mat, show.diag)
  } else if (type == "upper") {
    corr <- .get_upper_tri(corr, show.diag)
    p.mat <- .get_upper_tri(p.mat, show.diag)
  }
  
  # Melt corr and pmat
  corr <- reshape2::melt(corr, na.rm = TRUE, as.is = as.is)
  colnames(corr) <- c("Var1", "Var2", "value")
  corr$pvalue <- rep(NA, nrow(corr))
  corr$signif <- rep(NA, nrow(corr))
  
  if (!is.null(p.mat)) {
    p.mat <- reshape2::melt(p.mat, na.rm = TRUE)
    corr$coef <- corr$value
    corr$pvalue <- p.mat$value
    corr$signif <- as.numeric(p.mat$value <= sig.level)
    p.mat <- subset(p.mat, p.mat$value < sig.level)
    if (insig == "blank") {
      corr$value <- corr$value * corr$signif
    }
  }
  
  
  corr$abs_corr <- abs(corr$value) * 10
  
  # heatmap
  p <-
    ggplot2::ggplot(
      data = corr,
      mapping = ggplot2::aes_string(x = "Var1", y = "Var2", fill = "value")
    )
  
  # modification based on method
  if (method == "square") {
    p <- p +
      ggplot2::geom_tile(color = outline.color)
  } else if (method == "circle") {
    p <- p +
      ggplot2::geom_point(
        color = outline.color,
        shape = 21,
        ggplot2::aes_string(size = "abs_corr")
      ) +
      ggplot2::scale_size(range = c(4, 10)) +
      ggplot2::guides(size = "none")
  }
  
  # adding colors
  p <- p + ggplot2::scale_fill_gradient2(
    low = colors[1],
    high = colors[3],
    mid = colors[2],
    midpoint = 0,
    limit = c(-1, 1),
    space = "Lab",
    name = legend.title
  )
  
  # depending on the class of the object, add the specified theme
  if (class(ggtheme)[[1]] == "function") {
    p <- p + ggtheme()
  } else if (class(ggtheme)[[1]] == "theme") {
    p <- p + ggtheme
  }
  
  
  p <- p +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = tl.srt,
        vjust = 1,
        size = tl.cex,
        hjust = 1
      ),
      axis.text.y = ggplot2::element_text(size = tl.cex)
    ) +
    ggplot2::coord_fixed()
  
  label <- round(x = corr[, "value"], digits = digits)
  if (!is.null(p.mat) & insig == "blank") {
    ns <- corr$pvalue > sig.level
    if (sum(ns) > 0) label[ns] <- " "
  }
  
  # matrix cell labels
  if (lab) {
    p <- p +
      ggplot2::geom_text(
        mapping = ggplot2::aes_string(x = "Var1", y = "Var2"),
        label = label,
        color = lab_col,
        size = lab_size
      )
  }
  
  # matrix cell glyphs
  if (!is.null(p.mat) & insig == "pch") {
    p <- p + ggplot2::geom_point(
      data = p.mat,
      mapping = ggplot2::aes_string(x = "Var1", y = "Var2"),
      shape = pch,
      size = pch.cex,
      color = pch.col
    )
  }
  
  # add titles
  if (title != "") {
    p <- p +
      ggplot2::ggtitle(title)
  }
  
  # removing legend
  if (!show.legend) {
    p <- p +
      ggplot2::theme(legend.position = "none")
  }
  
  # removing panel
  p <- p +
    .no_panel()
  p
}



#' Compute the matrix of correlation p-values
#'
#' @param x numeric matrix or data frame
#' @param ... other arguments to be passed to the function cor.test.
#' @rdname ggcorrplot
#' @export

cor_pmat <- function(x, ...) {
  
  # initializing values
  mat <- as.matrix(x)
  n <- ncol(mat)
  p.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  
  # creating the p-value matrix
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- stats::cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  
  # name rows and columns of the p-value matrix
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  
  # return the final matrix
  p.mat
}



#+++++++++++++++++++++++
# Helper Functions
#+++++++++++++++++++++++

# Get lower triangle of the correlation matrix
.get_lower_tri <- function(cormat, show.diag = FALSE) {
  if (is.null(cormat)) {
    return(cormat)
  }
  cormat[upper.tri(cormat)] <- NA
  if (!show.diag) {
    diag(cormat) <- NA
  }
  return(cormat)
}

# Get upper triangle of the correlation matrix
.get_upper_tri <- function(cormat, show.diag = FALSE) {
  if (is.null(cormat)) {
    return(cormat)
  }
  cormat[lower.tri(cormat)] <- NA
  if (!show.diag) {
    diag(cormat) <- NA
  }
  return(cormat)
}

.remove_diag <- function(cormat) {
  if (is.null(cormat)) {
    return(cormat)
  }
  diag(cormat) <- NA
  cormat
}
# hc.order correlation matrix
.hc_cormat_order <- function(cormat, hc.method = "complete") {
  dd <- stats::as.dist((1 - cormat) / 2)
  hc <- stats::hclust(dd, method = hc.method)
  hc$order
}

.no_panel <- function() {
  ggplot2::theme(
    axis.title.x = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank()
  )
}


# Convert a tbl to matrix
.tibble_to_matrix <- function(x) {
  x <- as.data.frame(x)
  rownames(x) <- x[, 1]
  x <- x[, -1]
  as.matrix(x)
}

#after cutting down the APM and csde1 outliers there are 176 samples, well above the 119 samples required for power.

# write all of the tmm data for these samples into a csv for further survival analysis

APMCSDE1outlierscut <-  as.data.frame(TCGAOCCAMSmergedComBatBatchC[,colnames(TCGAOCCAMSmergedComBatBatchC) %in% rownames(regcut6) ])

write.csv(APMCSDE1outlierscut, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/OCCAMSTCGATMMMergedBatchCorrectedAPMCSDE1OutliersCut.csv")



#boxplot genes of interest
TMMforBox <- read.csv( file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/OCCAMSTCGATMMMergedBatchCorrectedAPMCSDE1OutliersCut.csv", header = TRUE, row.names = 1)
library(tidyverse)
library(hrbrthemes)
library(viridis)

TMMforBox2 <- as.data.frame(t(TMMforBox))
TMMforBox3 <- as.data.frame(TMMforBox2$CSDE1)
TMMforBox3$Gene <- "CSDE1"
colnames(TMMforBox3) <- c("Gene_expression_TMM", "Gene")
TMMforBox3 %>%
  ggplot( aes(x=Gene, y=value, fill=Gene)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=1, alpha=0.9) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("A boxplot with jitter") +
  xlab("")

TMMforBox3 %>%
  ggplot( aes(x=Gene, y=Gene_expression_TMM, fill=Gene)) +
  geom_violin() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("CSDE1 gene expression") +
  xlab("")  +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=1, alpha=0.9) +
  theme_ipsum() 

#write out quantiles for later analysis if needed

quantilesCSDE1 <- quantile(TMMforBox2$CSDE1)
UpperquantileCSDE1 <- TMMforBox2[TMMforBox2$CSDE1 >518.13062,]
LowerquantileCSDE1 <- TMMforBox2[TMMforBox2$CSDE1 <317.50950,]
rownames(UpperquantileCSDE1) <- strtrim(rownames(UpperquantileCSDE1), 12)
rownames(LowerquantileCSDE1) <- strtrim(rownames(LowerquantileCSDE1), 12)

write.csv(UpperquantileCSDE1, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/UpperquantileCSDE1.csv")
write.csv(LowerquantileCSDE1, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/LowerquantileCSDE1.csv")

quantilesHLA.E <- quantile(TMMforBox2$'HLA-E')
UpperquantileHLA.E <- TMMforBox2[TMMforBox2$'HLA-E' >751.52543,]
LowerquantileHLA.E <- TMMforBox2[TMMforBox2$'HLA-E' <401.30496,]
rownames(UpperquantileHLA.E) <- strtrim(rownames(UpperquantileHLA.E), 12)
rownames(LowerquantileHLA.E) <- strtrim(rownames(LowerquantileHLA.E), 12)

write.csv(UpperquantileHLA.E, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/UpperquantileHLA.E.csv")
write.csv(LowerquantileHLA.E, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/LowerquantileHLA.E.csv")

quantilesTAPBP <- quantile(TMMforBox2$'TAPBP')
UpperquantileTAPBP <- TMMforBox2[TMMforBox2$'TAPBP' >392.14891,]
LowerquantileTAPBP <- TMMforBox2[TMMforBox2$'TAPBP' <237.45954,]
rownames(UpperquantileTAPBP) <- strtrim(rownames(UpperquantileTAPBP), 12)
rownames(LowerquantileTAPBP) <- strtrim(rownames(LowerquantileTAPBP), 12)

write.csv(UpperquantileTAPBP, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/UpperquantileTAPBP.csv")
write.csv(LowerquantileTAPBP, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/LowerquantileTAPBP.csv")

quantilesPSMB9 <- quantile(TMMforBox2$'PSMB9')
UpperquantilePSMB9 <- TMMforBox2[TMMforBox2$'PSMB9' >53.155879 ,]
LowerquantilePSMB9 <- TMMforBox2[TMMforBox2$'PSMB9' <22.644066  ,]
rownames(UpperquantilePSMB9) <- strtrim(rownames(UpperquantilePSMB9), 12)
rownames(LowerquantilePSMB9) <- strtrim(rownames(LowerquantilePSMB9), 12)

write.csv(UpperquantilePSMB9, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/UpperquantilePSMB9.csv")
write.csv(LowerquantilePSMB9, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/LowerquantilePSMB9.csv")

quantilesSPPL2A <- quantile(TMMforBox2$'SPPL2A')
UpperquantileSPPL2A <- TMMforBox2[TMMforBox2$'SPPL2A' >87.88697 ,]
LowerquantileSPPL2A <- TMMforBox2[TMMforBox2$'SPPL2A' <52.72214  ,]
rownames(UpperquantileSPPL2A) <- strtrim(rownames(UpperquantileSPPL2A), 12)
rownames(LowerquantileSPPL2A) <- strtrim(rownames(LowerquantileSPPL2A), 12)

write.csv(UpperquantileSPPL2A, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/UpperquantileSPPL2A.csv")
write.csv(LowerquantileSPPL2A, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/LowerquantileSPPL2A.csv")

quantilesHLA.DRB1 <- quantile(TMMforBox2$'HLA-DRB1')
UpperquantileHLA.DRB1 <- TMMforBox2[TMMforBox2$'HLA-DRB1' >182.46198  ,]
LowerquantileHLA.DRB1 <- TMMforBox2[TMMforBox2$'HLA-DRB1' <59.48811  ,]
rownames(UpperquantileHLA.DRB1) <- strtrim(rownames(UpperquantileHLA.DRB1), 12)
rownames(LowerquantileHLA.DRB1) <- strtrim(rownames(LowerquantileHLA.DRB1), 12)

write.csv(UpperquantileHLA.DRB1, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/UpperquantileHLA.DRB1.csv")
write.csv(LowerquantileHLA.DRB1, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/LowerquantileHLA.DRB1.csv")

quantilesHLA.DQA1 <- quantile(TMMforBox2$'HLA-DQA1')
UpperquantileHLA.DQA1 <- TMMforBox2[TMMforBox2$'HLA-DQA1' >41.43338   ,]
LowerquantileHLA.DQA1 <- TMMforBox2[TMMforBox2$'HLA-DQA1' <13.19299    ,]
rownames(UpperquantileHLA.DQA1) <- strtrim(rownames(UpperquantileHLA.DQA1), 12)
rownames(LowerquantileHLA.DQA1) <- strtrim(rownames(LowerquantileHLA.DQA1), 12)

write.csv(UpperquantileHLA.DQA1, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/UpperquantileHLA.DQA1.csv")
write.csv(LowerquantileHLA.DQA1, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/LowerquantileHLA.DQA1.csv")

quantilesCIITA <- quantile(TMMforBox2$'CIITA')
UpperquantileCIITA <- TMMforBox2[TMMforBox2$'CIITA' >74.972055    ,]
LowerquantileCIITA <- TMMforBox2[TMMforBox2$'CIITA' <24.063123      ,]
rownames(UpperquantileCIITA) <- strtrim(rownames(UpperquantileCIITA), 12)
rownames(LowerquantileCIITA) <- strtrim(rownames(LowerquantileCIITA), 12)

write.csv(UpperquantileCIITA, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/UpperquantileCIITA.csv")
write.csv(LowerquantileCIITA, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/LowerquantileCIITA.csv")

quantilesMR1 <- quantile(TMMforBox2$'MR1')
UpperquantileMR1 <- TMMforBox2[TMMforBox2$'MR1' >25.087139     ,]
LowerquantileMR1 <- TMMforBox2[TMMforBox2$'MR1' <12.764532      ,]
rownames(UpperquantileMR1) <- strtrim(rownames(UpperquantileMR1), 12)
rownames(LowerquantileMR1) <- strtrim(rownames(LowerquantileMR1), 12)

write.csv(UpperquantileMR1, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/UpperquantileMR1.csv")
write.csv(LowerquantileMR1, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/LowerquantileMR1.csv")

quantilesCD1D <- quantile(TMMforBox2$'CD1D')
UpperquantileCD1D <- TMMforBox2[TMMforBox2$'CD1D' >1.9615490    ,]
LowerquantileCD1D <- TMMforBox2[TMMforBox2$'CD1D' <0.7472583          ,]
rownames(UpperquantileCD1D) <- strtrim(rownames(UpperquantileCD1D), 12)
rownames(LowerquantileCD1D) <- strtrim(rownames(LowerquantileCD1D), 12)

write.csv(UpperquantileCD1D, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/UpperquantileCD1D.csv")
write.csv(LowerquantileCD1D, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/LowerquantileCD1D.csv")

quantilesERAP2 <- quantile(TMMforBox2$'ERAP2')
UpperquantileERAP2 <- TMMforBox2[TMMforBox2$'ERAP2' >93.583610     ,]
LowerquantileERAP2 <- TMMforBox2[TMMforBox2$'ERAP2' <17.935874           ,]
rownames(UpperquantileERAP2) <- strtrim(rownames(UpperquantileERAP2), 12)
rownames(LowerquantileERAP2) <- strtrim(rownames(LowerquantileERAP2), 12)

write.csv(UpperquantileERAP2, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/UpperquantileERAP2.csv")
write.csv(LowerquantileERAP2, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/LowerquantileERAP2.csv")

quantilesCALR <- quantile(TMMforBox2$'CALR')
UpperquantileCALR <- TMMforBox2[TMMforBox2$'CALR' >638.3340       ,]
LowerquantileCALR <- TMMforBox2[TMMforBox2$'CALR' <388.2636            ,]
rownames(UpperquantileCALR) <- strtrim(rownames(UpperquantileCALR), 12)
rownames(LowerquantileCALR) <- strtrim(rownames(LowerquantileCALR), 12)

write.csv(UpperquantileCALR, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/UpperquantileCALR.csv")
write.csv(LowerquantileCALR, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/LowerquantileCALR.csv")

quantilesPSMB8 <- quantile(TMMforBox2$'PSMB8')
UpperquantilePSMB8 <- TMMforBox2[TMMforBox2$'PSMB8' >108.864826   ,]
LowerquantilePSMB8 <- TMMforBox2[TMMforBox2$'PSMB8' <54.854514       ,]
rownames(UpperquantilePSMB8) <- strtrim(rownames(UpperquantilePSMB8), 12)
rownames(LowerquantilePSMB8) <- strtrim(rownames(LowerquantilePSMB8), 12)

write.csv(UpperquantilePSMB8, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/UpperquantilePSMB8.csv")
write.csv(LowerquantilePSMB8, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/LowerquantilePSMB8.csv")

quantilesPSMB9 <- quantile(TMMforBox2$'PSMB9')
UpperquantilePSMB9 <- TMMforBox2[TMMforBox2$'PSMB9' >53.155879    ,]
LowerquantilePSMB9 <- TMMforBox2[TMMforBox2$'PSMB9' < 22.644066     ,]
rownames(UpperquantilePSMB9) <- strtrim(rownames(UpperquantilePSMB9), 12)
rownames(LowerquantilePSMB9) <- strtrim(rownames(LowerquantilePSMB9), 12)

write.csv(UpperquantilePSMB9, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/UpperquantilePSMB9.csv")
write.csv(LowerquantilePSMB9, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/LowerquantilePSMB9.csv")

quantilesHLA.DRA <- quantile(TMMforBox2$'HLA-DRA')
UpperquantileHLA.DRA <- TMMforBox2[TMMforBox2$'HLA-DRA' >430.27959   ,]
LowerquantileHLA.DRA <- TMMforBox2[TMMforBox2$'HLA-DRA' <147.18064       ,]
rownames(UpperquantileHLA.DRA) <- strtrim(rownames(UpperquantileHLA.DRA), 12)
rownames(LowerquantileHLA.DRA) <- strtrim(rownames(LowerquantileHLA.DRA), 12)

write.csv(UpperquantileHLA.DRA, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/UpperquantileHLA.DRA.csv")
write.csv(LowerquantileHLA.DRA, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/LowerquantileHLA.DRA.csv")

quantilesHLA.DPA1 <- quantile(TMMforBox2$'HLA-DPA1')
UpperquantileHLA.DPA1 <- TMMforBox2[TMMforBox2$'HLA-DPA1' >166.53841    ,]
LowerquantileHLA.DPA1 <- TMMforBox2[TMMforBox2$'HLA-DPA1' <51.71401   ,]
rownames(UpperquantileHLA.DPA1) <- strtrim(rownames(UpperquantileHLA.DPA1), 12)
rownames(LowerquantileHLA.DPA1) <- strtrim(rownames(LowerquantileHLA.DPA1), 12)

write.csv(UpperquantileHLA.DPA1, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/UpperquantileHLA.DPA1.csv")
write.csv(LowerquantileHLA.DPA1, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/LowerquantileHLA.DPA1.csv")

quantilesLGMN <- quantile(TMMforBox2$'LGMN')
UpperquantileLGMN <- TMMforBox2[TMMforBox2$'LGMN' >117.533527 ,]
LowerquantileLGMN <- TMMforBox2[TMMforBox2$'LGMN' <67.031446  ,]
rownames(UpperquantileLGMN) <- strtrim(rownames(UpperquantileLGMN), 12)
rownames(LowerquantileLGMN) <- strtrim(rownames(LowerquantileLGMN), 12)

write.csv(UpperquantileLGMN, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/UpperquantileLGMN.csv")
write.csv(LowerquantileLGMN, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/LowerquantileLGMN.csv")

quantilesCTSL <- quantile(TMMforBox2$'CTSL')
UpperquantileCTSL <- TMMforBox2[TMMforBox2$'CTSL' >68.953180       ,]
LowerquantileCTSL <- TMMforBox2[TMMforBox2$'CTSL' <25.800936      ,]
rownames(UpperquantileCTSL) <- strtrim(rownames(UpperquantileCTSL), 12)
rownames(LowerquantileCTSL) <- strtrim(rownames(LowerquantileCTSL), 12)

write.csv(UpperquantileCTSL, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/UpperquantileCTSL.csv")
write.csv(LowerquantileCTSL, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/LowerquantileCTSL.csv")

quantilesPSMB10 <- quantile(TMMforBox2$'PSMB10')
UpperquantilePSMB10 <- TMMforBox2[TMMforBox2$'PSMB10' >16.5640787        ,]
LowerquantilePSMB10 <- TMMforBox2[TMMforBox2$'PSMB10' <8.7328417       ,]
rownames(UpperquantilePSMB10) <- strtrim(rownames(UpperquantilePSMB10), 12)
rownames(LowerquantilePSMB10) <- strtrim(rownames(LowerquantilePSMB10), 12)

write.csv(UpperquantilePSMB10, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/UpperquantilePSMB10.csv")
write.csv(LowerquantilePSMB10, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/LowerquantilePSMB10.csv")

quantilesCD74 <- quantile(TMMforBox2$'CD74')
UpperquantileCD74 <- TMMforBox2[TMMforBox2$'CD74' >1091.8056         ,]
LowerquantileCD74 <- TMMforBox2[TMMforBox2$'CD74' <430.8701         ,]
rownames(UpperquantileCD74) <- strtrim(rownames(UpperquantileCD74), 12)
rownames(LowerquantileCD74) <- strtrim(rownames(LowerquantileCD74), 12)

write.csv(UpperquantileCD74, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/UpperquantileCD74.csv")
write.csv(LowerquantileCD74, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/LowerquantileCD74.csv")


quantilesHLA.DQA1 <- quantile(TMMforBox2$'HLA-DQA1')
UpperquantileHLA.DQA1 <- TMMforBox2[TMMforBox2$'HLA-DQA1' >41.43338          ,]
LowerquantileHLA.DQA1 <- TMMforBox2[TMMforBox2$'HLA-DQA1' <13.19299           ,]
rownames(UpperquantileHLA.DQA1) <- strtrim(rownames(UpperquantileHLA.DQA1), 12)
rownames(LowerquantileHLA.DQA1) <- strtrim(rownames(LowerquantileHLA.DQA1), 12)

write.csv(UpperquantileHLA.DQA1, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/UpperquantileHLA.DQA1.csv")
write.csv(LowerquantileHLA.DQA1, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/LowerquantileHLA.DQA1.csv")


