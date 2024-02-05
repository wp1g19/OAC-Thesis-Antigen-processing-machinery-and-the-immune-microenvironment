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
dge <- calcNormFactors(dge,method = "TMM", lib.size=NULL)
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
y <- calcNormFactors(y, method = "TMM")

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
pca <- pca(t(OCCAMSTMMMatrixnoNA), ncomp = 3)


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


pca <- pca(t(OCCAMSTMMMatrixnoNAoutliercut), ncomp = 3)
colorlist <- rainbow(17)
Scatter_Density(pca, batch = batch, color.set = colorlist)
#After removing poor quality samples no batch effect between OCCAMS sites was discovered.

#correlate expression of CSDE1 and APM genes in the occams data set

cor1 <- as.matrix(t(OCCAMSTMMMatrixnoNAoutliercut))

csdgenes <- c("HLA-A", "HLA-B", "HLA-C", "HLA-E", "TAP1", "TAP2", "TAPBP", "B2M", "CSDE1")

cor1 <- OCCAMSTMMMatrixnoNAoutliercutBatch[rownames(OCCAMSTMMMatrixnoNAoutliercutBatch) %in% csdgenes, ]
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
regOC <- as.data.frame(t(OCCAMSTMMMatrixnoNAoutliercutBatch))

quantile(regOC$`HLA-A`, probs = seq(.1, .9, by = .1))
quantile(regOC$CSDE1, probs = seq(.1, .9, by = .1))
regOC1 <- regOC[regOC$CSDE1 >50,]
regOC2 <- regOC1[regOC1$CSDE1 <1000,]
regOC3 <- regOC2[regOC2$`HLA-A` <3000,]
regOC4 <- regOC3[regOC3$`HLA-A` >300,]

ggplot(regOC4, aes(x=CSDE1, y=`HLA-A`)) + 
  geom_point() +stat_smooth(method = "lm", col = "red") +
  stat_cor(label.x = 3, label.y = 2000) 



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





