library(TCGAbiolinks)
library(recount)
library(biomaRt)
library(CePa)
library(dplyr)
# Download and preprocess GTex data

########Query from Recount2 platform: esophageal Carcinoma#######
#ESCA.recount.gtex<-TCGAquery_recount2(project="GTEX", tissue="esophagus")
ESCA.recount.gtex <- read.gct(file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/GTEx/gene_reads_stomach.gct")
ESCA.recount.gtex <- as.data.frame(ESCA.recount.gtex)
ESCA.recount.gtex <- ESCA.recount.gtex[!duplicated(ESCA.recount.gtex$Description),]
rownames(ESCA.recount.gtex) <- ESCA.recount.gtex$Description
ESCA.recount.gtex$Description <- NULL
ESCA.recount.gtex <- mutate_all(ESCA.recount.gtex, function(x) as.numeric(as.character(x)))
ESCA.recount.gtex <-as.matrix(ESCA.recount.gtex)

library(edgeR)
dge <- DGEList(counts = ESCA.recount.gtex)
keep  <- filterByExpr(dge)
table(keep)
dge <- dge[keep, , keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge,method = "TMM", lib.size=NULL)
GTexTMMMatrix <- cpm(dge)
GTexTMMMatrix <- as.data.frame(GTexTMMMatrix)

#read in the batch corrected tcga/occams data
tcgaoccamsTMM <- read.csv(file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/OCCAMSTCGATMMMergedBatchCorrectedAPMCSDE1OutliersCut.csv", header = TRUE, row.names = 1)

tcgaoccamsTMMfilt <- tcgaoccamsTMM[row.names(tcgaoccamsTMM) %in% row.names(GTexTMMMatrix), ]
GTexTMMMatrixfilt <- GTexTMMMatrix[row.names(GTexTMMMatrix) %in% row.names(tcgaoccamsTMMfilt), ]

tcgaocaamGtexmergedTMM <- merge(tcgaoccamsTMMfilt, GTexTMMMatrixfilt, by = 0)
rownames(tcgaocaamGtexmergedTMM) <- tcgaocaamGtexmergedTMM$Row.names
tcgaocaamGtexmergedTMM$Row.names <- NULL

# batch correction
library(pcaExplorer)
library(mixOmics)
library(pcaExplorer)
library(PLSDAbatch)
batch <- c(rep("TCGAOCCAMS",176),
           rep("GTex", 359)
           )

pca <- mixOmics::pca(t(tcgaocaamGtexmergedTMM), ncomp = 3)
colorlist <- rainbow(2)
Scatter_Density(pca, batch = batch, color.set = colorlist)
library(sva)

tcgaocaamGtexmergedTMMCombatBatch <- ComBat(
  tcgaocaamGtexmergedTMM,
  batch ,
  mod = NULL,
  par.prior = TRUE,
  prior.plots = FALSE,
  mean.only = FALSE,
  ref.batch = NULL,
  BPPARAM = bpparam("SerialParam")
)
pca <- mixOmics::pca(t(tcgaocaamGtexmergedTMMCombatBatch), ncomp = 3)
colorlist <- rainbow(2)
Scatter_Density(pca, batch = batch, color.set = colorlist)

#split out the dataframes into the tumour and normal for DEA
tcgaocaamGtexmergedTMMCombatBatch[tcgaocaamGtexmergedTMMCombatBatch < 0] <- 0
TCGAOCCAMSsplit <- tcgaocaamGtexmergedTMMCombatBatch[,1:176]
GTexsplit <- tcgaocaamGtexmergedTMMCombatBatch[,177:375]
tcgaocaamGtexmergedTMMCombatBatch <- as.data.frame(tcgaocaamGtexmergedTMMCombatBatch)
write.csv(tcgaocaamGtexmergedTMMCombatBatch, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/BatchedTCGAOCCAMSGTEX.csv")
# Diff.expr.analysis (DEA)

sampledata <- read.csv(file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/Sample data OCCAMTCGAGTEX.csv", header = TRUE)

# Create the DEseq2DataSet object
library(DESeq2)
tcgaocaamGtexmergedTMMCombatBatch <- mutate_all(tcgaocaamGtexmergedTMMCombatBatch, function(x) as.integer(as.numeric(x)))
sampledata$Sample.type <- factor(sampledata$Sample.type, levels=c("Primary_Tumour", "Normal"))
deseq2Data <- DESeqDataSetFromMatrix(countData=tcgaocaamGtexmergedTMMCombatBatch, colData=sampledata, design= ~Sample.type)
dim(deseq2Data)
dim(deseq2Data[rowSums(counts(deseq2Data)) > 5, ])
# Perform pre-filtering of the data
deseq2Data <- deseq2Data[rowSums(counts(deseq2Data)) > 0, ]

#if (!requireNamespace("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")
#BiocManager::install("BiocParallel", version = "3.8")

# Register the number of cores to use
library(BiocParallel)
register(SnowParam(4))
# 1. Run pipeline for differential expression steps (if you set up parallel processing, set parallel = TRUE here)
deseq2Data <- DESeq(deseq2Data, parallel = TRUE)

# Extract differential expression results
# For "tissueType" perform primary vs normal comparison
deseq2Results <- results(deseq2Data, contrast= c("Sample.type", "Primary_Tumour", "Normal"))

# View summary of results
summary(deseq2Results)                         

# Using DEseq2 built in method
plotMA(deseq2Results)

# Generate logical column 
deseq2_res_all <- data.frame(deseq2Results) %>% mutate(threshold = padj < 0)

# Create the volcano plot
library(ggplot2)
ggplot2::ggplot(deseq2_res_all) + 
  ggplot2::geom_point(ggplot2::aes(x = log2FoldChange, y = -log10(padj), color = threshold)) + 
  ggplot2::xlab("log2 fold change") + 
  ggplot2::ylab("-log10 adjusted p-value") + 
  ggplot2::theme(legend.position = "none", 
        plot.title = ggplot2::element_text(size = ggplot2::rel(1.5), hjust = 0.5), 
        axis.title = ggplot2::element_text(size = ggplot2::rel(1.25)))
#if (!requireNamespace('BiocManager', quietly = TRUE))
#  install.packages('BiocManager')

#BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)
EnhancedVolcano:: EnhancedVolcano(deseq2_res_all,
                                  lab = rownames(deseq2_res_all),
                                  x = 'log2FoldChange',
                                  y = 'pvalue',
                                  title = 'TCGA/OCCAMS vs GTex Normal Stomach',
                                  pCutoff = 10e-6,
                                  FCcutoff = 0.5,
                                  pointSize = 3.0,
                                  labSize = 6.0)

#table out the genes of interest
genes <- c("NLRC5",
           "CSDE1",
           "CIITA",
           "CANX",
           "ERAP1",
           "CD1C",
           "HLA-DRB1",
           "ERAP2",
           "CD1A",
           "SPPL2A",
           "TAPBP",
           "CD1B",
           "HLA-DQA1",
           "PDIA3",
           "LGMN",
           "CTSS",
           "CD1D",
           "PSMB9",
           "MR1",
           "TAPBPL",
           "TAP2",
           "HLA-A",
           "HLA-DQA2",
           "HLA-B",
           "HLA-DPA1",
           "TAP1",
           "HLA-DQB2",
           "HLA-DMA",
           "HLA-DRB5",
           "CALR",
           "HLA-DPB1",
           "B2M",
           "CD74",
           "RFXANK",
           "RFX5",
           "HLA-DRA",
           "HLA-C",
           "HLA-G",
           "PSMB8",
           "HLA-DOA",
           "IRF1",
           "CTSL",
           "HLA-DMB",
           "HLA-E",
           "PSMB10",
           "RFXAP",
           "NEK2")
deseq2_res_all$Gene <- rownames(deseq2_res_all)
apmDEA <- deseq2_res_all[deseq2_res_all$Gene %in% genes ,]


#TCGA alone
library("stringr") 
#read in the batch corrected tcga/occams data
tcgaoccamsTMM <- read.csv(file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/OCCAMSTCGATMMMergedBatchCorrectedAPMCSDE1OutliersCut.csv", header = TRUE, row.names = 1)

tcgaoccamsTMMfilt <- tcgaoccamsTMM[row.names(tcgaoccamsTMM) %in% row.names(GTexTMMMatrix), ]
GTexTMMMatrixfilt <- GTexTMMMatrix[row.names(GTexTMMMatrix) %in% row.names(tcgaoccamsTMMfilt), ]
tcgaoccamsTMMfilt <- tcgaoccamsTMMfilt[,str_detect(colnames(tcgaoccamsTMMfilt), "TCGA") ]
tcgaocaamGtexmergedTMM <- merge(tcgaoccamsTMMfilt, GTexTMMMatrixfilt, by = 0)
rownames(tcgaocaamGtexmergedTMM) <- tcgaocaamGtexmergedTMM$Row.names
tcgaocaamGtexmergedTMM$Row.names <- NULL

# batch correction
library(pcaExplorer)
library(mixOmics)
library(pcaExplorer)
library(PLSDAbatch)
batch <- c(rep("TCGAOCCAMS",69),
           rep("GTex", 375)
           )

pca <- mixOmics::pca(t(tcgaocaamGtexmergedTMM), ncomp = 3)
colorlist <- rainbow(2)
Scatter_Density(pca, batch = batch, color.set = colorlist)
library(sva)

tcgaocaamGtexmergedTMMCombatBatch <- ComBat(
  tcgaocaamGtexmergedTMM,
  batch ,
  mod = NULL,
  par.prior = TRUE,
  prior.plots = FALSE,
  mean.only = FALSE,
  ref.batch = NULL,
  BPPARAM = bpparam("SerialParam")
)
pca <- mixOmics::pca(t(tcgaocaamGtexmergedTMMCombatBatch), ncomp = 3)
colorlist <- rainbow(2)
Scatter_Density(pca, batch = batch, color.set = colorlist)

#split out the dataframes into the tumour and normal for DEA
tcgaocaamGtexmergedTMMCombatBatch[tcgaocaamGtexmergedTMMCombatBatch < 0] <- 0
TCGAOCCAMSsplit <- tcgaocaamGtexmergedTMMCombatBatch[,1:69]
GTexsplit <- tcgaocaamGtexmergedTMMCombatBatch[,69:375]
tcgaocaamGtexmergedTMMCombatBatch <- as.data.frame(tcgaocaamGtexmergedTMMCombatBatch)
write.csv(tcgaocaamGtexmergedTMMCombatBatch, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/BatchedTCGAGTEX.csv")
# Diff.expr.analysis (DEA)

sampledata <- read.csv(file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/Sample data TCGAGTEX.csv", header = TRUE)

# Create the DEseq2DataSet object
library(DESeq2)
library(dplyr)
tcgaocaamGtexmergedTMMCombatBatch <- mutate_all(tcgaocaamGtexmergedTMMCombatBatch, function(x) as.integer(as.numeric(x)))
sampledata$Sample.type <- factor(sampledata$Sample.type, levels=c("Primary_Tumour", "Normal"))
deseq2Data <- DESeqDataSetFromMatrix(countData=tcgaocaamGtexmergedTMMCombatBatch, colData=sampledata, design= ~Sample.type)
dim(deseq2Data)
dim(deseq2Data[rowSums(counts(deseq2Data)) > 5, ])
# Perform pre-filtering of the data
deseq2Data <- deseq2Data[rowSums(counts(deseq2Data)) > 0, ]

#if (!requireNamespace("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")
#BiocManager::install("BiocParallel", version = "3.8")

# Register the number of cores to use
library(BiocParallel)
register(SnowParam(4))
# 1. Run pipeline for differential expression steps (if you set up parallel processing, set parallel = TRUE here)
deseq2Data <- DESeq(deseq2Data, parallel = TRUE)

# Extract differential expression results
# For "tissueType" perform primary vs normal comparison
deseq2Results <- results(deseq2Data, contrast= c("Sample.type", "Primary_Tumour", "Normal"))

# View summary of results
summary(deseq2Results)                         

# Using DEseq2 built in method
plotMA(deseq2Results)

# Generate logical column 
deseq2_res_all <- data.frame(deseq2Results) %>% mutate(threshold = padj < 1)

# Create the volcano plot
library(ggplot2)
ggplot2::ggplot(deseq2_res_all) + 
  ggplot2::geom_point(ggplot2::aes(x = log2FoldChange, y = -log10(padj), color = threshold)) + 
  ggplot2::xlab("log2 fold change") + 
  ggplot2::ylab("-log10 adjusted p-value") + 
  ggplot2::theme(legend.position = "none", 
        plot.title = ggplot2::element_text(size = ggplot2::rel(1.5), hjust = 0.5), 
        axis.title = ggplot2::element_text(size = ggplot2::rel(1.25)))
#if (!requireNamespace('BiocManager', quietly = TRUE))
#  install.packages('BiocManager')

#BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)
EnhancedVolcano:: EnhancedVolcano(deseq2_res_all,
                                  lab = rownames(deseq2_res_all),
                                  x = 'log2FoldChange',
                                  y = 'pvalue',
                                  title = 'TCGA vs GTex Normal Gastro-oesophageal junction',
                                  pCutoff = 10e-6,
                                  FCcutoff = 0.5,
                                  pointSize = 3.0,
                                  labSize = 6.0)


#table out the genes of interest
genes <- c("NLRC5",
           "CSDE1",
           "CIITA",
           "CANX",
           "ERAP1",
           "CD1C",
           "HLA-DRB1",
           "ERAP2",
           "CD1A",
           "SPPL2A",
           "TAPBP",
           "CD1B",
           "HLA-DQA1",
           "PDIA3",
           "LGMN",
           "CTSS",
           "CD1D",
           "PSMB9",
           "MR1",
           "TAPBPL",
           "TAP2",
           "HLA-A",
           "HLA-DQA2",
           "HLA-B",
           "HLA-DPA1",
           "TAP1",
           "HLA-DQB2",
           "HLA-DMA",
           "HLA-DRB5",
           "CALR",
           "HLA-DPB1",
           "B2M",
           "CD74",
           "RFXANK",
           "RFX5",
           "HLA-DRA",
           "HLA-C",
           "HLA-G",
           "PSMB8",
           "HLA-DOA",
           "IRF1",
           "CTSL",
           "HLA-DMB",
           "HLA-E",
           "PSMB10",
           "RFXAP")
deseq2_res_all$Gene <- rownames(deseq2_res_all)
TCGAapmDEA <- deseq2_res_all[deseq2_res_all$Gene %in% genes ,]


#OCCAMS alone
library("stringr") 
#read in the batch corrected tcga/occams data
tcgaoccamsTMM <- read.csv(file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/OCCAMSTCGATMMMergedBatchCorrectedAPMCSDE1OutliersCut.csv", header = TRUE, row.names = 1)

tcgaoccamsTMMfilt <- tcgaoccamsTMM[row.names(tcgaoccamsTMM) %in% row.names(GTexTMMMatrix), ]
GTexTMMMatrixfilt <- GTexTMMMatrix[row.names(GTexTMMMatrix) %in% row.names(tcgaoccamsTMMfilt), ]
tcgaoccamsTMMfilt <- tcgaoccamsTMMfilt[,str_detect(colnames(tcgaoccamsTMMfilt), "OC.") ]
tcgaocaamGtexmergedTMM <- merge(tcgaoccamsTMMfilt, GTexTMMMatrixfilt, by = 0)
rownames(tcgaocaamGtexmergedTMM) <- tcgaocaamGtexmergedTMM$Row.names
tcgaocaamGtexmergedTMM$Row.names <- NULL

# batch correction
library(pcaExplorer)
library(mixOmics)
library(pcaExplorer)
library(PLSDAbatch)
batch <- c(rep("TCGAOCCAMS",107),
           rep("GTex", 375)
)

pca <- mixOmics::pca(t(tcgaocaamGtexmergedTMM), ncomp = 3)
colorlist <- rainbow(2)
Scatter_Density(pca, batch = batch, color.set = colorlist)
library(sva)

tcgaocaamGtexmergedTMMCombatBatch <- ComBat(
  tcgaocaamGtexmergedTMM,
  batch ,
  mod = NULL,
  par.prior = TRUE,
  prior.plots = FALSE,
  mean.only = FALSE,
  ref.batch = NULL,
  BPPARAM = bpparam("SerialParam")
)
pca <- mixOmics::pca(t(tcgaocaamGtexmergedTMMCombatBatch), ncomp = 3)
colorlist <- rainbow(2)
Scatter_Density(pca, batch = batch, color.set = colorlist)

#split out the dataframes into the tumour and normal for DEA
tcgaocaamGtexmergedTMMCombatBatch[tcgaocaamGtexmergedTMMCombatBatch < 0] <- 0
TCGAOCCAMSsplit <- tcgaocaamGtexmergedTMMCombatBatch[,1:69]
GTexsplit <- tcgaocaamGtexmergedTMMCombatBatch[,69:375]
tcgaocaamGtexmergedTMMCombatBatch <- as.data.frame(tcgaocaamGtexmergedTMMCombatBatch)
write.csv(tcgaocaamGtexmergedTMMCombatBatch, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/BatchedOCCAMSGTEX.csv")
# Diff.expr.analysis (DEA)

sampledata <- read.csv(file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/Sample data OCCAMGTEX.csv", header = TRUE)

# Create the DEseq2DataSet object
library(DESeq2)
library(dplyr)
tcgaocaamGtexmergedTMMCombatBatch <- mutate_all(tcgaocaamGtexmergedTMMCombatBatch, function(x) as.integer(as.numeric(x)))
sampledata$Sample.type <- factor(sampledata$Sample.type, levels=c("Primary_Tumour", "Normal"))
deseq2Data <- DESeqDataSetFromMatrix(countData=tcgaocaamGtexmergedTMMCombatBatch, colData=sampledata, design= ~Sample.type)
dim(deseq2Data)
dim(deseq2Data[rowSums(counts(deseq2Data)) > 5, ])
# Perform pre-filtering of the data
deseq2Data <- deseq2Data[rowSums(counts(deseq2Data)) > 0, ]

#if (!requireNamespace("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")
#BiocManager::install("BiocParallel", version = "3.8")

# Register the number of cores to use
library(BiocParallel)
register(SnowParam(4))
# 1. Run pipeline for differential expression steps (if you set up parallel processing, set parallel = TRUE here)
deseq2Data <- DESeq(deseq2Data, parallel = TRUE)

# Extract differential expression results
# For "tissueType" perform primary vs normal comparison
deseq2Results <- results(deseq2Data, contrast= c("Sample.type", "Primary_Tumour", "Normal"))

# View summary of results
summary(deseq2Results)                         

# Using DEseq2 built in method
plotMA(deseq2Results)

# Generate logical column 
deseq2_res_all <- data.frame(deseq2Results) %>% mutate(threshold = padj < 1)

# Create the volcano plot
library(ggplot2)
ggplot2::ggplot(deseq2_res_all) + 
  ggplot2::geom_point(ggplot2::aes(x = log2FoldChange, y = -log10(padj), color = threshold)) + 
  ggplot2::xlab("log2 fold change") + 
  ggplot2::ylab("-log10 adjusted p-value") + 
  ggplot2::theme(legend.position = "none", 
                 plot.title = ggplot2::element_text(size = ggplot2::rel(1.5), hjust = 0.5), 
                 axis.title = ggplot2::element_text(size = ggplot2::rel(1.25)))
#if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

#BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)
EnhancedVolcano:: EnhancedVolcano(deseq2_res_all,
                                  lab = rownames(deseq2_res_all),
                                  x = 'log2FoldChange',
                                  y = 'pvalue',
                                  title = 'OCCAMS vs GTex Normal Gastro-oesophageal junction',
                                  pCutoff = 10e-6,
                                  FCcutoff = 0.5,
                                  pointSize = 3.0,
                                  labSize = 6.0)

#table out the genes of interest
genes <- c("NLRC5",
           "CSDE1",
           "CIITA",
           "CANX",
           "ERAP1",
           "CD1C",
           "HLA-DRB1",
           "ERAP2",
           "CD1A",
           "SPPL2A",
           "TAPBP",
           "CD1B",
           "HLA-DQA1",
           "PDIA3",
           "LGMN",
           "CTSS",
           "CD1D",
           "PSMB9",
           "MR1",
           "TAPBPL",
           "TAP2",
           "HLA-A",
           "HLA-DQA2",
           "HLA-B",
           "HLA-DPA1",
           "TAP1",
           "HLA-DQB2",
           "HLA-DMA",
           "HLA-DRB5",
           "CALR",
           "HLA-DPB1",
           "B2M",
           "CD74",
           "RFXANK",
           "RFX5",
           "HLA-DRA",
           "HLA-C",
           "HLA-G",
           "PSMB8",
           "HLA-DOA",
           "IRF1",
           "CTSL",
           "HLA-DMB",
           "HLA-E",
           "PSMB10",
           "RFXAP")
deseq2_res_all$Gene <- rownames(deseq2_res_all)
OCCAMSapmDEA <- deseq2_res_all[deseq2_res_all$Gene %in% genes ,]


#OAC vs OSCC
esca.subtype_OAC <- TCGAquery_subtype("esca")
subset_OSCC_filter <- subset(esca.subtype_OAC, esca.subtype_OAC$`Histological Type - Oesophagus` == 'ESCC')
query.exp <- GDCquery(project = "TCGA-ESCA", 
                      legacy = FALSE,
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "STAR - Counts", sample.type = "Primary Tumor", barcode = subset_OSCC_filter$patient
                      )
GDCdownload(query.exp)
OSCCESCA.exp <- GDCprepare(query = query.exp, save = TRUE, save.filename = "OACESCAExp.rda")
library(dplyr)

dataPrep <- TCGAanalyze_Preprocessing(object = OSCCESCA.exp)


library(biomaRt)
dataPrep <- as.data.frame(dataPrep)
dataPrep$ensembl <- gsub("\\..*","", rownames(dataPrep))
dataPrep = dataPrep[!duplicated(dataPrep$ensembl),]
rownames(dataPrep) <- dataPrep$ensembl
dataPrep$ensembl <- NULL
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
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
TCGAOSCCTMMMatrix <- cpm(dge)
TCGAOSCCTMMMatrix <- as.data.frame(TCGAOSCCTMMMatrix)
TCGAOSCCTMMMatrix <- TCGAOSCCTMMMatrix[row.names(TCGAOSCCTMMMatrix) %in% row.names(tcgaoccamsTMM), ]
tcgaoccamsTMM <- tcgaoccamsTMM[row.names(tcgaoccamsTMM) %in% row.names(TCGAOSCCTMMMatrix), ]


mergedTCGAOCCAMSOAC_OSCC <- cbind(tcgaoccamsTMM, TCGAOSCCTMMMatrix)
pca <- mixOmics::pca(t(mergedTCGAOCCAMSOAC_OSCC), ncomp = 3)
colorlist <- rainbow(3)
batch <- c(rep("TCGAOAC",69),
           rep("OCCAMS", 107),
           rep("TCGAOSCC", 90)
)
Scatter_Density(pca, batch = batch, color.set = colorlist)
library(sva)


mergedTCGAOCCAMSOAC_OSCCTMMCombatBatch <- ComBat(
  mergedTCGAOCCAMSOAC_OSCC,
  batch ,
  mod = NULL,
  par.prior = TRUE,
  prior.plots = FALSE,
  mean.only = FALSE,
  ref.batch = NULL,
  BPPARAM = bpparam("SerialParam")
)
pca <- mixOmics::pca(t(mergedTCGAOCCAMSOAC_OSCCTMMCombatBatch), ncomp = 3)
colorlist <- rainbow(3)
Scatter_Density(pca, batch = batch, color.set = colorlist)

#split out the dataframes into the tumour and normal for DEA
mergedTCGAOCCAMSOAC_OSCCTMMCombatBatch[mergedTCGAOCCAMSOAC_OSCCTMMCombatBatch < 0] <- 0
TCGAOCCAMSsplit <- mergedTCGAOCCAMSOAC_OSCCTMMCombatBatch[,1:176]
TCGAOACsplit <- mergedTCGAOCCAMSOAC_OSCCTMMCombatBatch[,177:255]
tcgaocaamGtexmergedTMMCombatBatch <- as.data.frame(tcgaocaamGtexmergedTMMCombatBatch)
write.csv(mergedTCGAOCCAMSOAC_OSCCTMMCombatBatch, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/BatchedTCGAOCCAMS_OAC_OSCC.csv")
# Diff.expr.analysis (DEA)

sampledata <- read.csv(file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/Sample data OCCAMTCGA_OAC_OSCC.csv", header = TRUE)

# Create the DEseq2DataSet object
library(DESeq2)
library(dplyr)
mergedTCGAOCCAMSOAC_OSCCTMMCombatBatch <- as.data.frame(mergedTCGAOCCAMSOAC_OSCCTMMCombatBatch)
mergedTCGAOCCAMSOAC_OSCCTMMCombatBatch <- mutate_all(mergedTCGAOCCAMSOAC_OSCCTMMCombatBatch, function(x) as.integer(as.numeric(x)))
sampledata$Sample.type <- factor(sampledata$Sample.type, levels=c("OAC", "OSCC"))
deseq2Data <- DESeqDataSetFromMatrix(countData=mergedTCGAOCCAMSOAC_OSCCTMMCombatBatch, colData=sampledata, design= ~Sample.type)
dim(deseq2Data)
dim(deseq2Data[rowSums(counts(deseq2Data)) > 5, ])
# Perform pre-filtering of the data
deseq2Data <- deseq2Data[rowSums(counts(deseq2Data)) > 0, ]

#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("BiocParallel", version = "3.8")

# Register the number of cores to use
library(BiocParallel)
register(SnowParam(4))
# 1. Run pipeline for differential expression steps (if you set up parallel processing, set parallel = TRUE here)
deseq2Data <- DESeq(deseq2Data, parallel = TRUE)

# Extract differential expression results
# For "tissueType" perform primary vs normal comparison
deseq2Results <- results(deseq2Data, contrast= c("Sample.type", "OAC", "OSCC"))

# View summary of results
summary(deseq2Results)                         

# Using DEseq2 built in method
plotMA(deseq2Results)

# Generate logical column 
deseq2_res_all <- data.frame(deseq2Results) %>% mutate(threshold = padj < 1)

# Create the volcano plot
library(ggplot2)
ggplot2::ggplot(deseq2_res_all) + 
  ggplot2::geom_point(ggplot2::aes(x = log2FoldChange, y = -log10(padj), color = threshold)) + 
  ggplot2::xlab("log2 fold change") + 
  ggplot2::ylab("-log10 adjusted p-value") + 
  ggplot2::theme(legend.position = "none", 
                 plot.title = ggplot2::element_text(size = ggplot2::rel(1.5), hjust = 0.5), 
                 axis.title = ggplot2::element_text(size = ggplot2::rel(1.25)))
#if (!requireNamespace('BiocManager', quietly = TRUE))
#install.packages('BiocManager')

#BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)
EnhancedVolcano:: EnhancedVolcano(deseq2_res_all,
                                  lab = rownames(deseq2_res_all),
                                  x = 'log2FoldChange',
                                  y = 'pvalue',
                                  title = 'OAC vs OSCC',
                                  pCutoff = 0.05,
                                  FCcutoff = 0.5,
                                  pointSize = 3.0,
                                  labSize = 6.0)

#table out the genes of interest
genes <- c("NLRC5",
           "CSDE1",
           "CIITA",
           "CANX",
           "ERAP1",
           "CD1C",
           "HLA-DRB1",
           "ERAP2",
           "CD1A",
           "SPPL2A",
           "TAPBP",
           "CD1B",
           "HLA-DQA1",
           "PDIA3",
           "LGMN",
           "CTSS",
           "CD1D",
           "PSMB9",
           "MR1",
           "TAPBPL",
           "TAP2",
           "HLA-A",
           "HLA-DQA2",
           "HLA-B",
           "HLA-DPA1",
           "TAP1",
           "HLA-DQB2",
           "HLA-DMA",
           "HLA-DRB5",
           "CALR",
           "HLA-DPB1",
           "B2M",
           "CD74",
           "RFXANK",
           "RFX5",
           "HLA-DRA",
           "HLA-C",
           "HLA-G",
           "PSMB8",
           "HLA-DOA",
           "IRF1",
           "CTSL",
           "HLA-DMB",
           "HLA-E",
           "PSMB10",
           "RFXAP")
deseq2_res_all$Gene <- rownames(deseq2_res_all)
OCCAMSapmDEA <- deseq2_res_all[deseq2_res_all$Gene %in% genes ,]

