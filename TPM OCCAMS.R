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

genelengths <- as.data.frame(combined_df$genelength)
rownames(genelengths) <- rownames(combined_df)

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


#convert gene IDs
library(biomaRt)
combined_df_removedcols <- as.data.frame(genelengths)
combined_df_removedcols$gene <- rownames(combined_df_removedcols)
rownames(combined_df_removedcols) <- gsub("\\..*","", rownames(combined_df_removedcols))
#combined_df_removedcols = combined_df_removedcols[!duplicated(rownames(combined_df_removedcols)),]
mart <- useEnsembl(biomart = "ensembl", 
                   dataset = "hsapiens_gene_ensembl", 
                   mirror = "useast")
genes <- rownames(combined_df_removedcols)
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
combined_df_removedcols2 <- as.data.frame(combined_df_removedcols[ rownames(combined_df_removedcols) %in% G_list$ensembl_gene_id ,])
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
genelengths <- mutate_all(combined_df_removedcols, function(x) as.numeric(as.character(x)))
genelengths$gene <- NULL

#install.packages("remotes")
#remotes::install_github("davidrequena/drfun")
library(DRnaSeq)
OCCAMSTPM <- tpm(OCCAMSMatrixnoNA, gene_lengths = genelengths$`combined_df$genelength`)

#Read in TCGA TPM data
library(TCGAbiolinks)
library(SummarizedExperiment)
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
TCGARawcounts <- combined_df_removedcols
library(DRnaSeq)
combined_df_removedcolsfilt <- as.data.frame(combined_df_removedcols[rownames(OCCAMSTPM),])
TCGATPM <- tpm(combined_df_removedcolsfilt, gene_lengths = genelengths$`combined_df$genelength`)

colnames(TCGATPM) <- strtrim(colnames(TCGATPM), 12)

barcodes <- c("TCGA.L7.A6VZ", "TCGA.IG.A4QS", "TCGA.L5.A8NV", "TCGA.R6.A6XG", "TCGA.R6.A6DQ", "TCGA.L5.A8NE", "TCGA.L5.A8NT", "TCGA.L5.A4OQ", "TCGA.2H.A9GH", "TCGA.R6.A8W8", "TCGA.L5.A8NJ", "TCGA.V5.AASX", "TCGA.L5.A4OP", "TCGA.L5.A8NF", "TCGA.L5.A8NM", "TCGA.L5.A8NR", "TCGA.2H.A9GN", "TCGA.L5.A8NS", "TCGA.L5.A4OI", "TCGA.2H.A9GO", "TCGA.R6.A8W5", "TCGA.L5.A88V", "TCGA.L5.A8NI", "TCGA.L5.A4OW", "TCGA.L5.A4OJ", "TCGA.L5.A8NW", "TCGA.R6.A6L6", "TCGA.L5.A891", "TCGA.RE.A7BO", "TCGA.IC.A6RE", "TCGA.JY.A6F8", "TCGA.JY.A93E", "TCGA.VR.AA4D", "TCGA.R6.A6KZ", "TCGA.S8.A6BV", "TCGA.X8.AAAR", "TCGA.L5.A4OO", "TCGA.JY.A6FB", "TCGA.JY.A93C", "TCGA.2H.A9GK", "TCGA.2H.A9GF", "TCGA.2H.A9GI", "TCGA.2H.A9GL", "TCGA.2H.A9GM", "TCGA.2H.A9GR","TCGA.JY.A6FH", "TCGA.JY.A93D", "TCGA.JY.A938", "TCGA.L5.A4OE", "TCGA.L5.A4OF", "TCGA.L5.A4OG","TCGA.L5.A4OH", "TCGA.L5.A4ON", "TCGA.L5.A4OR", "TCGA.L5.A4OS", "TCGA.L5.A4OT", "TCGA.L5.A4OU", "TCGA.L5.A4OX", "TCGA.L5.A8NG", "TCGA.L5.A8NH", "TCGA.L5.A8NL", "TCGA.L5.A8NN", "TCGA.L5.A8NU", "TCGA.L5.A43C", "TCGA.L5.A43I", "TCGA.L5.A88Y", "TCGA.M9.A5M8", "TCGA.Q9.A6FW", "TCGA.R6.A6DN", "TCGA.R6.A6L4", "TCGA.R6.A6XQ", "TCGA.R6.A6Y2", "TCGA.R6.A8WC", "TCGA.R6.A8WG", "TCGA.V5.A7RB", "TCGA.V5.A7RE", "TCGA.V5.AASW", "TCGA.VR.A8EQ", "TCGA.ZR.A9CJ")
colnames(TCGATPM) <- gsub("-",".", colnames(TCGATPM))
TCGATPMfilt <- as.data.frame(TCGATPM[, colnames(TCGATPM) %in% barcodes])

TCGATPMfilt2 <- as.data.frame(TCGATPMfilt[rownames(TCGATPMfilt) %in% rownames(OCCAMSTPM), ])
OCCAMSTPMfilt <- as.data.frame(OCCAMSTPM[rownames(OCCAMSTPM) %in% rownames(TCGATPMfilt2), ])
#ID conversion for OCCAMS
# Convert the SLX to OCCAMS number
OCCAMSTPMMatrixnoNA2 <- t(OCCAMSTPMfilt)
conversionmap <- read.csv("C:/Users/wp1g19/Desktop/Post return writing/TPR/Cibersortx/SLXtoIDOCCAMS.csv")
rownames(conversionmap) <- conversionmap$SLX_ID
rownames(conversionmap) <- gsub("-",".", rownames(conversionmap))
rownames(OCCAMSTPMMatrixnoNA2) <- gsub("-",".", rownames(OCCAMSTPMMatrixnoNA2))
rownames(conversionmap) <- gsub("_",".", rownames(conversionmap))
rownames(OCCAMSTPMMatrixnoNA2) <- gsub("_",".", rownames(OCCAMSTPMMatrixnoNA2))

OCCAMSTPMMatrixnoNA3 <- as.data.frame(OCCAMSTPMMatrixnoNA2)
rownames(OCCAMSTPMMatrixnoNA3) = substr(rownames(OCCAMSTPMMatrixnoNA3),1,nchar(rownames(OCCAMSTPMMatrixnoNA3))-13)
OCCAMSTPMMatrixnoNA3 <- as.data.frame(OCCAMSTPMMatrixnoNA3)
OCCAMSTPMMatrixnoNA3 <- OCCAMSTPMMatrixnoNA3[ rownames(OCCAMSTPMMatrixnoNA3) %in% rownames(conversionmap),]
conversionmaporder <- conversionmap[row.names(OCCAMSTPMMatrixnoNA3) ,]

row.names(OCCAMSTPMMatrixnoNA3) <- conversionmaporder$ID
OCCAMSTPMfilt <- as.data.frame(t(OCCAMSTPMMatrixnoNA3))


TCGATPMfilt2$gene <- rownames(TCGATPMfilt2)
OCCAMSTPMfilt$gene <- rownames(OCCAMSTPMfilt)

OCCAMSTCGATPM <- merge(TCGATPMfilt2, OCCAMSTPMfilt, by = "gene")
rownames(OCCAMSTCGATPM) <- OCCAMSTCGATPM$gene
OCCAMSTCGATPM$gene <- NULL
write.csv(OCCAMSTCGATPM, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Digital cytometry/TPM rerun/TCGAOCCAMSTPM.csv")

batch <- c(rep("TCGA",79), rep("OCCAMS", 151))
library(pcaExplorer)
library(mixOmics)
library(pcaExplorer)
library(PLSDAbatch)
pca <- mixOmics::pca(t(OCCAMSTCGATPM), ncomp = 3)
Scatter_Density(pca)
colorlist <- rainbow(2)
Scatter_Density(pca, batch = batch, color.set = colorlist)

OCCAMSTCGATMMMatrixnoNAoutliercut <- OCCAMSTCGATPM
#[,!colnames(OCCAMSTCGATPM) %in% outlier_barcode ]

batch <- c(rep("TCGA",79), rep("OCCAMS", 151))
pca <- mixOmics::pca(t(OCCAMSTCGATMMMatrixnoNAoutliercut), ncomp = 3)
Scatter_Density(pca)
colorlist <- rainbow(2)
Scatter_Density(pca, batch = batch, color.set = colorlist)

library(sva)
batch <- c(rep("TCGA",79), rep("OCCAMS", 151))
OCCAMSTCGATPMCOMBAT <- sva::ComBat(OCCAMSTCGATPM,
                                            batch = batch,
                                            prior.plots = TRUE,
                                            mean.only = FALSE)
library(rrcov)
pcaHub <- PcaHubert(t(OCCAMSTCGATPMCOMBAT))
outliers <- which(pcaHub@flag=='FALSE')
outliers <- as.data.frame(pcaHub@flag)
filter <- split(outliers, ~outliers$`pcaHub@flag`)
keepcases <- filter[["TRUE"]]
OCCAMSTCGATPMCOMBATOUTLIERCUT <- OCCAMSTCGATPMCOMBAT[,colnames(OCCAMSTCGATPMCOMBAT) %in% rownames(keepcases) ]

write.csv(OCCAMSTCGATPMCOMBATOUTLIERCUT, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Digital cytometry/TPM rerun/TCGAOCCAMSTPM_BATCH_CORRECTED_OutliersCUt.csv")

pca <- mixOmics::pca(t(OCCAMSTCGATPMCOMBAT), ncomp = 3)
Scatter_Density(pca)
colorlist <- rainbow(2)
Scatter_Density(pca, batch = batch, color.set = colorlist)

library(edgeR)

OCCAMSTCGATPMEDGER <- removeBatchEffect(OCCAMSTCGATPM,
                                   batch = batch)
write.csv(OCCAMSTCGATPMCOMBAT, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Digital cytometry/TPM rerun/TCGAOCCAMSTPM_BATCH_CORRECTED.csv")
write.csv(OCCAMSTCGATPMEDGER, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Digital cytometry/TPM rerun/TCGAOCCAMSTPM_LIMMA_BATCH_CORRECTED.csv")

#Batch correct counts then tpm

# Convert the SLX to OCCAMS number
OCCAMSTPMMatrixnoNA2 <- t(OCCAMSMatrixnoNA)
conversionmap <- read.csv("C:/Users/wp1g19/Desktop/Post return writing/TPR/Cibersortx/SLXtoIDOCCAMS.csv")
rownames(conversionmap) <- conversionmap$SLX_ID
rownames(conversionmap) <- gsub("-",".", rownames(conversionmap))
rownames(OCCAMSTPMMatrixnoNA2) <- gsub("-",".", rownames(OCCAMSTPMMatrixnoNA2))
rownames(conversionmap) <- gsub("_",".", rownames(conversionmap))
rownames(OCCAMSTPMMatrixnoNA2) <- gsub("_",".", rownames(OCCAMSTPMMatrixnoNA2))

OCCAMSTPMMatrixnoNA3 <- as.data.frame(OCCAMSTPMMatrixnoNA2)
rownames(OCCAMSTPMMatrixnoNA3) = substr(rownames(OCCAMSTPMMatrixnoNA3),1,nchar(rownames(OCCAMSTPMMatrixnoNA3))-13)
OCCAMSTPMMatrixnoNA3 <- as.data.frame(OCCAMSTPMMatrixnoNA3)
OCCAMSTPMMatrixnoNA3 <- OCCAMSTPMMatrixnoNA3[ rownames(OCCAMSTPMMatrixnoNA3) %in% rownames(conversionmap),]
conversionmaporder <- conversionmap[row.names(OCCAMSTPMMatrixnoNA3) ,]

row.names(OCCAMSTPMMatrixnoNA3) <- conversionmaporder$ID
OCCAMSIDcorrected <- as.data.frame(t(OCCAMSTPMMatrixnoNA3))
colnames(TCGARawcounts) <- strtrim(colnames(TCGARawcounts), 12)
OCCAMSIDcorrected <- OCCAMSIDcorrected[ rownames(OCCAMSIDcorrected) %in% rownames(TCGARawcounts),]
TCGARawcounts <- TCGARawcounts[ rownames(TCGARawcounts) %in% rownames(OCCAMSIDcorrected),]
OCCAMSIDcorrected$gene <- rownames(OCCAMSIDcorrected)
TCGARawcounts$gene <- rownames(TCGARawcounts)

MergeRawCounts <- merge(OCCAMSIDcorrected, TCGARawcounts, by = "gene")
rownames(MergeRawCounts) <- MergeRawCounts$gene
MergeRawCounts$gene <- NULL


library(sva)
batch <- c(rep("OCCAMS", 151), rep("TCGA",79))
BatchCorrectedRaw <- sva::ComBat(MergeRawCounts,
                                   batch = batch,
                                   prior.plots = TRUE,
                                   mean.only = FALSE)
BatchCorrectedRaw <- removeBatchEffect(BatchCorrectedRaw,
                                        batch = batch)

pca <- mixOmics::pca(t(BatchCorrectedRaw), ncomp = 3)
Scatter_Density(pca)
colorlist <- rainbow(2)
Scatter_Density(pca, batch = batch, color.set = colorlist)
genelengths$gene <- rownames(genelengths)
genelengths2 <- as.data.frame(genelengths[rownames(BatchCorrectedRaw),])
BatchCorrectedTPM <- tpm(BatchCorrectedRaw, gene_lengths = genelengths2$`combined_df$genelength`)
write.csv(BatchCorrectedTPM, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Digital cytometry/TPM rerun/TCGAOCCAMSTPM_LIMMA_BATCH_CORRECTED_ALT.csv")
colnames(BatchCorrectedTPM) <- gsub("-",".", colnames(BatchCorrectedTPM))
colnames(BatchCorrectedTPM) <- gsub("/",".", colnames(BatchCorrectedTPM))
#filter to genes of interest
siggenes <- c("RFX5",
              "CTSS",
              "CD1D",
              "MR1",
              "SPPL2A",
              "CD74",
              "CIITA",
              "PSMB10",
              "LGMN",
              "CTSL",
              "TAPBPL",
              "CALR",
              "ERAP1",
              "ERAP2",
              "HLA-E",
              "HLA-DPA1",
              "PSMB9",
              "HLA_B",
              "HLA-DRA",
              "PSMB8",
              "HLA-DRB5",
              "HLA-DQA1",
              "HLA-DRB1",
              "TAPBP",
              "HLA-G",
              "HLA-A",
              "CSDE1")


quantileforbatch <- as.data.frame(t(BatchCorrectedTPM))
quantilesCSDE1 <- quantile(quantileforbatch$CSDE1)
UpperquantileCSDE1 <- quantileforbatch[quantileforbatch$CSDE1 >151.00116,]
LowerquantileCSDE1 <- quantileforbatch[quantileforbatch$CSDE1 <92.84667,]

write.csv(UpperquantileCSDE1, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/UpperquantileCSDE1.csv")
write.csv(LowerquantileCSDE1, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/LowerquantileCSDE1.csv")

